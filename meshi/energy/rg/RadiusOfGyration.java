/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy.rg;

import meshi.molecularElements.atoms.AtomList;
import meshi.molecularElements.atoms.Atom;
import meshi.util.Updateable;

/**
 *
 */
public class RadiusOfGyration implements Updateable {
    private AtomList atomList;

    public final AtomList atoms() {
        return atomList;
    }

    private double invNumberOfAtoms;
    private int size;

    public final int size() {
        return size;
    }

    public Atom atomAt(int index) {
        return atomList.get(index);
    }
    //
    private double[] dRGdX, dRGdY, dRGdZ;

    public final double dRGdX(int index) {
        return dRGdX[index];
    }

    public final double dRGdY(int index) {
        return dRGdY[index];
    }

    public final double dRGdZ(int index) {
        return dRGdZ[index];
    }
    //
    private CenterOfMass centerOfMass;
    private double dRGdCMx, dRGdCMy, dRGdCMz;
    private double invCMnumberOfAtoms;

    public final double dRGdCMX() {
        return dRGdCMx;
    }

    public final double dRGdCMy() {
        return dRGdCMy;
    }

    public final double dRGdCMz() {
        return dRGdCMz;
    }
    //
    private double radiusOfGyration, dRadiusOfGyration;

    public double radiusOfGyration() {
        return radiusOfGyration;
    }

    public double dRadiusOfGyration() {
        return dRadiusOfGyration;
    }
    //
    /*
     * Used to avoid reupdate in the same minimization step
     */
    private int numberOfUpdates = 0;


    public RadiusOfGyration(AtomList atomList, CenterOfMass centerOfMass) {
        this.centerOfMass = centerOfMass;
        invCMnumberOfAtoms = centerOfMass.invSize();
        //
        this.atomList = atomList;
        size = atomList.size();
        invNumberOfAtoms = 1.00 / size;
        //
        dRGdX = new double[size];
        dRGdY = new double[size];
        dRGdZ = new double[size];
    }

    public void update(int numberOfUpdates) {
        if (numberOfUpdates == this.numberOfUpdates + 1) {
            this.numberOfUpdates++;
            update();
        } else if (numberOfUpdates != this.numberOfUpdates)
            throw new RuntimeException("Something weird with RadiusOfGyration(int numberOfUpdates)\n" +
                    "numberOfUpdates = " + numberOfUpdates + " this.numberOfUpdates = " + this.numberOfUpdates);
    }

    private void update() {
        Atom atom;
        double sum, sigmaX = 0, sigmaY = 0, sigmaZ = 0;
        sum = 0;
        for (int i = 0; i < size; i++) {
            atom = atomList.get(i);
            dRGdX[i] = atom.x() - centerOfMass.getX();
            dRGdY[i] = atom.y() - centerOfMass.getY();
            dRGdZ[i] = atom.z() - centerOfMass.getZ();
            sum += dRGdX[i] * dRGdX[i] + dRGdY[i] * dRGdY[i] + dRGdZ[i] * dRGdZ[i];
            sigmaX += dRGdX[i];
            sigmaY += dRGdY[i];
            sigmaZ += dRGdZ[i];
        }
        radiusOfGyration = Math.sqrt(invNumberOfAtoms * sum);
        dRadiusOfGyration = 0.5 * invNumberOfAtoms / radiusOfGyration;
        for (int i = 0; i < size; i++) {
            dRGdX[i] = 2 * dRGdX[i] * dRadiusOfGyration;
            dRGdY[i] = 2 * dRGdY[i] * dRadiusOfGyration;
            dRGdZ[i] = 2 * dRGdZ[i] * dRadiusOfGyration;
        }
        dRGdCMx = invNumberOfAtoms * invCMnumberOfAtoms / radiusOfGyration * sigmaX;
        dRGdCMy = invNumberOfAtoms * invCMnumberOfAtoms / radiusOfGyration * sigmaY;
        dRGdCMz = invNumberOfAtoms * invCMnumberOfAtoms / radiusOfGyration * sigmaZ;
    }

    public void setNumberOfUpdates(int numberOfUpdates) {
        this.numberOfUpdates = numberOfUpdates;
    }

    public void test(Atom atom) {
       System.out.println(" RadiusOfGyration ");
       double delta = 0.00000001;
        update();

    }
}
