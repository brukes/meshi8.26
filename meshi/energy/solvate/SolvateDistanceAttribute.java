/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.solvate;

import meshi.util.*;

/**
 * In the SolvateEnergy evaluation function these fields are to be recalculated for EVERY
 * distance in the non-bonded list. They are used several times in the solvate evaluation
 * process, and we would like to calculate them once. To this end we attach this special
 * class as an attribute isOn each Distance instance.
 * For each distance we calculate two sigmoid values, and an hydrogen-bond strength value:
 * sigmCa1 - The carbon index of atomOne 2 isOn atomOne 1. If atomOne 2 is not a carbon then this value
 * should be zero. If atomOne 2 is a carbon then this value should be ~1.0 if atomOne 2
 * is spatially near atomOne 1. This index drops sigmoidally to zero the farther
 * atomOne 2 is.
 * sigmCa2 - The same as sigmCa1, except detailing the affect of atomOne 1 isOn atomOne 2.
 * sigmHBa1 - The hydrogen bond (HB) strength between atomOne 1 and 2. If atomOne 2 can not form HB
 * with atomOne 1, because of its chemical type then this value should be zero.
 * If atomOne 2 can form HB with atomOne 1 then this value should be ~1.0 if atomOne 1 and 2
 * are sufficiently close to create a HB, and if their orientation (defined
 * also by their base atoms - see HydrogenBondDahiyat) permit hydrogen bonding.
 * This value drops steeply to zero if the conditions to hydrogen bonding
 * are violated.
 * <p/>
 * Also provided are the carbon sigmoid values derivative relatives to the atomOne coordinates. They
 * have the general form: dsigmCa{1/2}d{x/y/z}{1/2}
 */


public class SolvateDistanceAttribute implements MeshiAttribute {


    public SolvateDistanceType type() {
        return null;
    }

    public int key() {
        return SOLVATE_ALL_ATOM_ATTRIBUTE;
    }

}











