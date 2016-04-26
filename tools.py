#!/usr/bin/env python2
# -*- coding: utf8 -*-
from os import listdir, path
from os.path import isfile, join
import sys
import glob
import fileinput
import re

list_src_files = ["module_param.f90", "module_frontvirt.f90",
                  "incompact3d.f90", "mesure.f90", "schemas.f90", "derive.f90",
                  "pression.f90", "spectral.f90", "outils.f90", "filtre.f90",
                  "parametre.f90", "front_virt.f90", "navier.f90",
                  "pression6.f90", "jmccm1d2.f90", "jmccm1d3.f90",
                  "jmccm1d4.f90", "jmccm1d5.f90", "jmccm1d.f90",
                  "jmccm1dp.f90", "jmcctranspcs.f90", "jmcsm1dxy.f90",
                  "jmerreur1.f90", "jmerreur2.f90", "jmfact.f90",
                  "jmfftfax.f90", "jmgetmessage.f90", "jmgetsetcode.f90",
                  "jmgetsetstop.f90", "jmgetstop.f90", "jmrfftmlt.f90",
                  "jmscm1dxy.f90", "jmsetcode.f90", "jmtable.f90",
                  "cosfftmlt.f90", "jmtableCos.f90", "jmtableSin.f90",
                  "slfft2d.f90"]


def look_for_implicit(fich):
    implicit_is_missing = False
    for line in open(fich):

        if re.search("^ *subroutine", line) is not None:
            if implicit_is_missing:
                print("missing implicit none if file " + fich +
                      " in subroutine " + fx)
            fx = line.split("subroutine ")[1].split("(")[0]
            implicit_is_missing = True

        if re.search("^ *function", line) is not None:
            if implicit_is_missing:
                print("missing implicit none if file " + fich +
                      " in subroutine " + fx)
            fx = line.split("function ")[1].split("(")[0]
            implicit_is_missing = True

        if "implicit none" in line:
            implicit_is_missing = False


def look_for_module(fich):
    for line in open(fich):
        if re.search("^ *module", line) is not None:
            return True
    return False


def look_for_use(fich):
    list = []
    for line in open(fich):
        if re.search("^ *call", line) is not None:
            elem = line.split("call")[1].split("(")[0].replace(" ", "")
            if elem not in list:
                list.append(elem)
    return list


def fx_to_file(fx):
    for fich in glob.glob('*.f90'):
        for line in open(fich):
            if re.search("^ *function", line) is not None:
                if fx in line:
                    return fich
            if re.search("^ *subroutine", line) is not None:
                if fx in line:
                    return fich
    return None


def add_uses(fich, list_use):
    for line in fileinput.input(fich, inplace=1):
        if fileinput.isfirstline():
            print(line),
            for use in list_use:
                print("  use " + use)
        else:
            print(line),


def add_module(fich):
    filename = path.splitext(fich)[0]
    module_name = filename+"_m"
    for line in fileinput.input(fich, inplace=1):
        if fileinput.isfirstline():
            print("module " + module_name)
            print("implicit none")
            print("contains")
        print(line),
    with open(fich, "a") as myfile:
        myfile.write("end module " + module_name)
    print("module " + module_name + " created in file " + fich)


for fich in glob.glob('*.f90'):
    look_for_implicit(fich)
    if not look_for_module(fich):
        add_module(fich)
    list_fx = look_for_use(fich)
    if list_fx:
        list_file = []
        for fx in list_fx:
            fich1 = fx_to_file(fx)
            if fich1 and not fich1 == fich and fich1 not in list_file:
                list_file.append(fich1)
        if list_file:
            list_use = []
            for fich1 in list_file:
                list_use.append(fich1.replace(".f90", "_m"))
            add_uses(fich, list_use)
#            print(fich),
#            print(list_use)



