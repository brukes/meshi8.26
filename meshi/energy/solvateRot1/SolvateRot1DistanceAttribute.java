/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvateRot1;

import meshi.util.MeshiAttribute;

/**
 * In the SolvateRot1Energy evaluation function these fields are to be recalculated for EVERY
 * distance in the non-bonded list. They are used several times in the solvate evaluation
 * process, and we would like to calculate them once. To this end we attach this special
 * class as an attribute isOn each Distance instance.
 * For each distance we calculate two sigmoid values: one for atomOne 1 in the distance (a1),
 * and one for atomOne 2 in the distance (a2):
 * sigmCa1 - The carbon index of atomOne 2 isOn atomOne 1. This value should be ~1.0 if atomOne 2
 * is spatially near atomOne 1. This index drops sigmoidally to zero the farther
 * atomOne 2 is.
 * sigmCa2 - The same as sigmCa1, except detailing the affect of atomOne 1 isOn atomOne 2.
 * <p/>
 * Also provided are the sigmoid values derivative relatives to the atomOne coordinates. They
 * have the general form: dsigmCa{1/2}d{x/y/z}{1/2}
 */


public class SolvateRot1DistanceAttribute implements MeshiAttribute {

    public SolvateRot1DistanceAttribute() {
    }

    public int key() {
        return SOLVATE_ROT1_ATTRIBUTE;
    }

    public double sigmCa1;
    public double sigmCa2;
    public double dsigmCa1dx1;
    public double dsigmCa1dy1;
    public double dsigmCa1dz1;
    public double dsigmCa1dx2;
    public double dsigmCa1dy2;
    public double dsigmCa1dz2;
    public double dsigmCa2dx1;
    public double dsigmCa2dy1;
    public double dsigmCa2dz1;
    public double dsigmCa2dx2;
    public double dsigmCa2dy2;
    public double dsigmCa2dz2;


    public final void resetAllSigmVals() {
        sigmCa1 = 0.0;
        dsigmCa1dx1 = 0.0;
        dsigmCa1dy1 = 0.0;
        dsigmCa1dz1 = 0.0;
        dsigmCa1dx2 = 0.0;
        dsigmCa1dy2 = 0.0;
        dsigmCa1dz2 = 0.0;
        sigmCa2 = 0.0;
        dsigmCa2dx1 = 0.0;
        dsigmCa2dy1 = 0.0;
        dsigmCa2dz1 = 0.0;
        dsigmCa2dx2 = 0.0;
        dsigmCa2dy2 = 0.0;
        dsigmCa2dz2 = 0.0;
    }

}