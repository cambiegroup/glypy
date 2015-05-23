
import operator

from pygly2.algorithms.database import RecordDatabase
from pygly2.utils import groupby
from .spectra.spectrum_model import MSMSSqlDB
from .spectra.decon2ls import parse_to_database, Decon2LSPeak
from .matching import MassShift, NoShift, make_struct, ppm_error


match_key_getter = operator.attrgetter("match_key")

PrecursorMatch = make_struct("PrecursorMatch", ("match_key", "mass", "ppm_error",
                                                "intensity", "charge", "spectrum_data", "scan_id"))


ResultsGroup = make_struct("ResultsGroup", ("match_key", "average_mass", "average_ppm_error",
                                            "total_intensity", "scan_count", "scan_density",
                                            "charge_state_count", "members"))


def density(x):
    try:
        minx, maxx = min(x), max(x)
        return len(x) / float(maxx - minx)
    except (ZeroDivisionError, ValueError):
        return 0


def combine_results(group):
    n = len(group)
    mass = 0
    ppm_error = 0
    intensity = 0
    charge_states = set()
    scans = set()
    for match in group:
        mass += match.mass
        ppm_error += match.ppm_error
        intensity += match.intensity
        charge_states.add(match.charge)
        scans.add(match.scan_id)
    average_mass = mass / float(n)
    average_ppm_error = ppm_error / float(n)
    charge_state_count = len(charge_states)
    scan_density = density(scans)
    grouped = ResultsGroup(
        group[0].match_key, average_mass, average_ppm_error,
        intensity, len(scans), scan_density, charge_state_count, group)
    return grouped


def match_decon2ls_isos(spectrum_db, hypothesis, adducts=None, ms1_match_tolerance=1e-5, output_path=None):
    if adducts is None:
        adducts = [NoShift]

    results = RecordDatabase(output_path, record_type=hypothesis.record_type)
    results.apply_schema()

    for composition in hypothesis:
        matches = []
        scans_searched = set()
        charge_states = set()
        for adduct in adducts:
            for spectrum in spectrum_db.ppm_match_tolerance_search(
                    composition.intact_mass + adduct.mass, ms1_match_tolerance):
                match_ppm = ppm_error(composition.intact_mass + adduct.mass, spectrum.neutral_mass)
                match = PrecursorMatch(
                    str(composition.id) + ":" + adduct.name, spectrum.neutral_mass, match_ppm,
                    spectrum.get("abundance"), spectrum.charge, spectrum.other_data,
                    spectrum.scan_ids[0])
                matches.append(match)
                scans_searched.update(spectrum.scan_ids)
                charge_states.add(spectrum.charge)

        groups = groupby(matches, match_key_getter)

        composition.precursor_matches = {key: combine_results(group) for key, group in groups.items()}
        composition.precursor_scans_searched = scans_searched
        composition.precursor_scan_density = density(scans_searched)
        composition.precursor_charge_states = charge_states
        results.load_data([composition], set_id=False, cast=True, commit=False)
    results.commit()
    results.apply_indices()
    return results
