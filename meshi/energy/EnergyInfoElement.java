/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.energy;

import meshi.util.info.DoubleInfoElement;
import meshi.util.info.DoubleInfoElementList;
import meshi.util.info.InfoType;
import meshi.util.info.MeshiInfoElement;

/**

 */
public class EnergyInfoElement extends DoubleInfoElement {
    private DoubleInfoElementList infoList;
    private double weight = Double.MIN_VALUE;
    public double weight() {
        return weight;
    }

    public EnergyInfoElement(InfoType type, String comment) {
        this(type, comment, 0,new DoubleInfoElementList());
    }

    public EnergyInfoElement(InfoType type, String comment, double weight) {
        this(type, comment, weight,new DoubleInfoElementList());
    }

    public EnergyInfoElement(InfoType type, String comment, double weight, DoubleInfoElementList infoList) {
        this(type, comment, weight,infoList,0);
    }
    public EnergyInfoElement(InfoType type,
                             String comment,
                             double weight,
                             DoubleInfoElementList infoList,
                             double value) {
        super(type, comment, value);
        this.weight = weight;
        this.infoList = infoList;
    }

    public EnergyInfoElement(EnergyInfoElement other) {
        this(other.type, other.comment, other.weight(),new DoubleInfoElementList(),other.value());
        if (other.infoList()!= null)  {
              for(DoubleInfoElement element : other.infoList())
                  infoList().add(element);
        }
    }

    public double energy() {
        return value();
    }

    public DoubleInfoElementList infoList() {
        return infoList;
    }

    public void resetEnergy() {
        setValue(0);
        if (infoList() != null) {
            for (MeshiInfoElement element : infoList()) {
                ((DoubleInfoElement) element).setValue(0);
            }
        }
    }
}
