# coding=windows-1251

__author__="artur"
__date__ ="$Dec 16, 2009 11:40:59 PM$"

from MechanicalFrame import *
from examples.ballandbeam.ballandbeam import bb_main
from examples.pendulums.simple_pendulum import sp_main
from examples.pendulums.gyro_pendulum import gp_main

def main():

    bb_main()
    #bb_template()
    #sp_main()
    #gp_main()

if __name__ == "__main__":
    main()
