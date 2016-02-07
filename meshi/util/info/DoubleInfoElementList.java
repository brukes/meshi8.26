/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.info;

import meshi.energy.EnergyInfoElement;

import java.util.ArrayList;

/**

 */
public class DoubleInfoElementList extends ArrayList<DoubleInfoElement> {
    DoubleInfoElement get(InfoType type) {
        DoubleInfoElement out = null;
        for (DoubleInfoElement element : this) {
            if (element.type.equals(type)) {
                if (out == null) out = element;
                else throw new RuntimeException("At least two elements of type " + type + " in list " + this);
            }
        }
        return out;
    }

    public DoubleInfoElement addElement(DoubleInfoElement newElement) {
        add(newElement);
        return newElement;
    }
}
