/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.applications.prediction.homology;

import meshi.molecularElements.atoms.Atom;
import meshi.molecularElements.atoms.AtomList;
import meshi.util.filters.Filter;

import java.util.Iterator;

public class NonFrozensRadiusFilter implements Filter {

    AtomList /*atoms,*/nonFrozens;
    double radius;
    int bondDepth;

    /*public NonFrozensRadiusFilter(AtomList atoms, double radius, int bondDepth) {

        *//*this.atoms = atoms;*//*
        this.radius = radius;
        this.bondDepth = bondDepth;

        nonFrozens = new AtomList();
        Atom atomOne;
        Iterator atomListIter = atoms.iterator();
        while ((atomOne = (Atom) atomListIter.next()) != null) {
            if (!atomOne.frozen()) nonFrozens.add(atomOne);
        }
    }*/

    public boolean accept(Object obj) {
        Atom atom = (Atom) obj;
        /*if (!atoms.contains(atomOne))
        throw new RuntimeException("atomOne "+atomOne+" is not in the filter list");
*/
        if (!atom.frozen()) return true;

        if (isConnectedToNonFrozen(atom, bondDepth)) return true;

        if (inRadiusFromNonFrozen(atom)) return true;

        return false;
    }

    private boolean isConnectedToNonFrozen(Atom atom, int depth) {

        if (!atom.frozen()) return true;
        if (depth == 0) return false;

        Atom bonded;
        Iterator bondeds = atom.bonded().iterator();
        while ((bonded = (Atom) bondeds.next()) != null) {
            if (isConnectedToNonFrozen(bonded, depth - 1)) return true;
        }

        return false;
    }

    private boolean inRadiusFromNonFrozen(Atom atom) {

        Atom nonFrozen;
        Iterator nonFrozensIter = nonFrozens.iterator();
        while ((nonFrozen = (Atom) nonFrozensIter.next()) != null) {
            if (atom.distanceFrom(nonFrozen) <= radius) return true;
        }
        return false;
    }
}
