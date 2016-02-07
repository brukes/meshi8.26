/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.pairwiseNonBondedTerms.CooperativeSumma;

import meshi.util.*;

public class SummaAttribute implements MeshiAttribute {

    public double halfEnergy, fx, fy, fz;

    public int key() {
        return SUMMA_ATTRIBUTE;
    }

}
