#!/usr/bin/env python3
#
# Copyright 2022-2023 ImpactX contributors
# Authors: Marco Garten
# License: BSD-3-Clause-LBNL
#
# -*- coding: utf-8 -*-

import warnings

from impactx import RefPart, elements

from .MADXParser import MADXParser


class sc:
    """
    This class is used in lieu of scipy.constants
    to avoid a direct dependency on it.
    At the time of writing, this file was the only
    one requiring scipy in the ImpactX source outside
    of examples.
    """

    c = 299792458.0
    electron_volt = 1.602176634e-19
    physical_constants = {"electron-muon mass ratio": (0.00483633169, "", 1.1e-10)}
    m_e = 9.1093837015e-31
    m_p = 1.67262192369e-27
    m_u = 1.6605390666e-27


def lattice(parsed_beamline, nslice=1):
    """
    Function that converts a list of elements in the MADXParser format into ImpactX format
    :param parsed_beamline: list of dictionaries
    :param nslice: number of ds slices per element
    :return: list of translated dictionaries
    """

    # mapping is "MAD-X": "ImpactX",
    madx_to_impactx_dict = {
        "MARKER": "None",
        "DRIFT": "Drift",
        "SBEND": "Sbend",  # Sector Bending Magnet
        "SOLENOID": "Sol",  # Ideal, thick Solenoid: MAD-X user guide 10.9 p78
        "QUADRUPOLE": "Quad",  # Quadrupole
        "DIPEDGE": "DipEdge",
        # Kicker, idealized thin element,
        # MADX also defines length "L" and a roll angle around the longitudinal axis "TILT"
        # https://mad.web.cern.ch/mad/webguide/manual.html#Ch11.S11
        "KICKER": "Kicker",  # TKICKER, HKICKER and VKICKER become KICKER elements
        # note: in MAD-X, this keeps track only of the beam centroid,
        # "In addition it serves to record the beam position for closed orbit correction."
        "MONITOR": "BeamMonitor",  # drift + output diagnostics
        "MULTIPOLE": "Multipole",
        "NLLENS": "NonlinearLens",
        "RFCAVITY": "RFCavity",
        # TODO Figure out how to identify these
        "ShortRF": "ShortRF",
        "ConstF": "ConstF",
    }

    impactx_beamline = []

    for d in parsed_beamline:
        # print(d)
        if d["type"] == "drift":
            impactx_beamline.append(elements.Drift(ds=d["l"], nslice=nslice))
        elif d["type"] == "quadrupole":
            if d["k1"] != 0.0:
                impactx_beamline.append(
                    elements.Quad(ds=d["l"], k=d["k1"], nslice=nslice)
                )
            else:
                impactx_beamline.append(
                    elements.Drift(ds=d["l"])
                )
                print("Warning: Zero strength quad not supported. Replaced with drift")
        elif d["type"] == "sbend":
            impactx_beamline.append(
                elements.Sbend(ds=d["l"], rc=d["l"] / d["angle"], nslice=nslice)
            )
        elif d["type"] == "solenoid":
            if d["ks"] != 0.0:
                impactx_beamline.append(
                    elements.Sol(ds=d["l"], ks=d["ks"], nslice=nslice)
                )
            else:
                impactx_beamline.append(
                    elements.Drift(ds=d["l"])
                )
                print("Warning: Zero strength solenoid not supported. Replaced with drift")
        elif d["type"] == "dipedge":
            impactx_beamline.append(
                elements.DipEdge(
                    psi=d["e1"],
                    rc=1.0 / d["h"],
                    # MAD-X is using half the gap height
                    g=2.0 * d["hgap"],
                    K2=d["fint"],
                )
            )
        elif d["type"] == "kicker":
            impactx_beamline.append(
                elements.Kicker(
                    xkick=d["hkick"],
                    ykick=d["vkick"],
                )
            )
        elif d["type"] == "vkicker":
            impactx_beamline.append(
                elements.Kicker(
                    xkick=0.0,
                    ykick=d['kick'],
                )
            )
        elif d["type"] == "hkicker":
            impactx_beamline.append(
                elements.Kicker(
                    xkick=d['kick'],
                    ykick=0.0,
                )
            )
        elif d["type"] == "monitor":
            if d["l"] > 0:
                impactx_beamline.append(elements.Drift(ds=d["l"], nslice=nslice))
            impactx_beamline.append(elements.BeamMonitor("monitor", backend="h5"))
        elif d["type"] == 'sextupole':
            print("Warning: Sextupole converted to drift")
            impactx_beamline.append(
                elements.Drift(ds=d['l'], nslice=nslice)
            )
        elif d["type"] == "marker":
            print("Warning: Marker element dropped")
            # TODO: marker listed in madx_to_impactx_dict but it looks like it is not supported here
        elif d["type"] == "rfcavity":
            impactx_beamline.append(
                #elements.RFCavity(
                #    ds=d["l"],
                #    escale=0.0,  # TODO: This cannot be set without reference particle information
                #    freq=d.get('freq', 0.0),  # If madx rfcavity uses harm this must be set after beamline created and with ref particle info
                #    phase=d['lag'],  # TODO: no idea if madx and impactx use same phase convention
                #    cos_coefficients=[0.0,],  # TODO: docs suggest these are optional and have defaults. Code disagrees.
                #    sin_coefficients=[0.0,]   # "
                #)
                elements.Drift(ds=d["l"])
            )
            if d.get('harmon'):
                print("Warning: Cavity freq cannot be set from harm automatically")

        else:
            raise NotImplementedError(
                "The beamline element named ",
                d["name"],
                "of type ",
                d["type"],
                "is not implemented in impactx.elements.",
                "Available elements are:",
                list(madx_to_impactx_dict.keys()),
            )
    return impactx_beamline

def generate_item_lookup(model, category):  # map _id to element dict
    assert category == 'elements' or category == 'beamlines'
    if category == 'elements':
        return {
            val['_id']: val for val in model.models.elements
        }
    elif category == 'beamlines':
        return {
            val['id']: val for val in model.models.beamlines
        }

def iterate_lattice(lattice_dict, beamline_map):
    for element in lattice_dict['items']:
        if element in beamline_map:
            new_beamline = beamline_map[element]
            for sub_element in iterate_lattice(new_beamline, beamline_map):
                yield sub_element
        else:
            # yield element _id
            yield element

def get_lattice_sirepo(madx_file, beamline=None):
    from sirepo import lib
    model = lib.Importer('madx').parse_file(madx_file)
    element_lookup = generate_item_lookup(model , 'elements')
    beamline_lookup = generate_item_lookup(model , 'beamlines')

    if beamline:
        for line in model.models.beamlines:
            if line.name == beamline:
                selected_line = line
                break
        else:
            raise NameError(f'requested beamline "{beamline}" not found in {madx_file}')
    else:
        selected_line = model.models.beamlines[0]

    flattened_line = iterate_lattice(selected_line, beamline_lookup)

    return [
        {**element_lookup[_id], 'type': element_lookup[_id]['type'].casefold()} for _id in flattened_line
    ]

def read_lattice(madx_file, nslice=1, beamline=None):
    """
    Function that reads elements from a MAD-X file into a list of ImpactX.KnownElements
    :param madx_file: file name to MAD-X file with beamline elements
    :param nslice: number of ds slices per element
    :return: list of ImpactX.KnownElements
    """
    # madx = MADXParser()
    #madx.parse(madx_file)
    # beamline = madx.getBeamline()
    # print(madx)
    beamline = get_lattice_sirepo(madx_file,beamline=beamline)
    return lattice(beamline, nslice)


def beam(particle, charge=None, mass=None, energy=None):
    """
    Function that converts a list of beam parameter dictionaries in the MADXParser format into ImpactX format

    Rules following https://mad.web.cern.ch/mad/releases/5.02.08/madxuguide.pdf pages 55f.

    :param str particle: reference particle name
    :param float charge: particle charge (proton charge units)
    :param float mass: particle mass (electron masses)
    :param float energy: total particle energy (GeV)
        - MAD-X default: 1 GeV
    :return dict: dictionary containing particle and beam attributes in ImpactX units
    """

    GeV2MeV = 1.0e3
    kg2MeV = sc.c**2 / sc.electron_volt * 1.0e-6
    muon_mass = sc.physical_constants["electron-muon mass ratio"][0] / sc.m_e
    if energy is None:
        energy_MeV = 1.0e3  # MAD-X default is 1 GeV total particle energy
    else:
        energy_MeV = energy * GeV2MeV

    impactx_beam = {
        "positron": {"mass": sc.m_e * kg2MeV, "charge": 1.0},
        "electron": {"mass": sc.m_e * kg2MeV, "charge": -1.0},
        "proton": {"mass": sc.m_p * kg2MeV, "charge": 1.0},
        "antiproton": {"mass": sc.m_p * kg2MeV, "charge": -1.0},
        "posmuon": {
            "mass": muon_mass * kg2MeV,
            "charge": 1.0,
        },  # positively charged muon (anti-muon)
        "negmuon": {
            "mass": muon_mass * kg2MeV,
            "charge": -1.0,
        },  # negatively charged muon
        "ion": {"mass": sc.m_u * kg2MeV, "charge": 1.0},
        "generic": {"mass": mass, "charge": charge},
    }

    if particle not in impactx_beam.keys():
        warnings.warn(
            f"Particle species name '{particle}' not in '{impactx_beam.keys()}'",
            UserWarning,
        )
        print(
            "Choosing generic particle species, using provided `charge`, `mass` and `energy`."
        )
        _particle = "generic"
    else:
        _particle = particle
        print(
            f"Choosing provided particle species '{particle}', ignoring potentially provided `charge` and `mass` and setting defaults.",
        )

    reference_particle = impactx_beam[_particle]
    reference_particle["energy"] = energy_MeV

    return reference_particle


def read_beam(ref: RefPart, madx_file):
    """
    Function that reads elements from a MAD-X file into a list of ImpactX.KnownElements
    :param RefPart ref: ImpactX reference particle (passed by reference)
    :param madx_file: file name to MAD-X file with beamline elements
    :return: list of ImpactX.KnownElements
    """
    # madx = MADXParser()
    # madx.parse(madx_file)
    from sirepo import lib
    model = lib.Importer('madx').parse_file(madx_file)
    for cmd in model.models.commands:
        if cmd._type == 'beam':
            madx_beam = cmd
            break
    else:
        raise Exception(f'No beam defined in MAD-X file {madx_file}')

    ref_particle_dict = beam(
        madx_beam['particle'],  # if particle species is known, mass, charge, and potentially energy are set to default
        # TODO MADX parser needs to extract charge if it's given,
        # TODO MADX parser needs to extract mass if it's given,
        energy=madx_beam['energy'],  # MADX default energy is in GeV
    )

    ref.set_charge_qe(ref_particle_dict["charge"])
    ref.set_mass_MeV(ref_particle_dict["mass"])
    ref.set_kin_energy_MeV(ref_particle_dict["energy"] - ref_particle_dict["mass"])

    return ref
