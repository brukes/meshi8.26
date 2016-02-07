/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.rg.filters;

import meshi.util.filters.Filter;
import meshi.molecularElements.atoms.Atom;
import meshi.parameters.SecondaryStructure;

/**
 *
 */
public class CoilFilter implements Filter {
    public boolean accept(Object o) {
        Atom atom = (Atom) o;
        if (atom.nowhere()) return false;
        SecondaryStructure secondaryStructure = atom.residue().secondaryStructure();
        if ((secondaryStructure == SecondaryStructure.HELIX) || (secondaryStructure == SecondaryStructure.SHEET))
            return false;
        if (secondaryStructure == SecondaryStructure.COIL) return true;
        throw new RuntimeException("Does not know what to do with atomOne " + atom + "\n" + "With residue secondary strctuure" + secondaryStructure);
    }

}
