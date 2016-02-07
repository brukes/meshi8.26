/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.info;

import java.util.ArrayList;

/**

 */
public class MeshiInfoElementList extends ArrayList<MeshiInfoElement> {

    public static MeshiInfoElementList deepCopy(MeshiInfoElementList list) {
        MeshiInfoElementList out = new MeshiInfoElementList();
        for (MeshiInfoElement element : list) {
            switch (element.valueType()) {
                case INTEGER:
                    out.add(new IntInfoElement((IntInfoElement) element));
                    break;
                case DOUBLE:
                    out.add(new DoubleInfoElement((DoubleInfoElement) element));
                    break;
                case STRING:
                    out.add(new StringInfoElement((StringInfoElement) element));
                    break;
                default:
                    throw new RuntimeException("Do not know how to handle " + element + "of value type" + element.valueType());
            }
        }
        return out;
    }

    public void add(MeshiInfoElementList list) {
        for (MeshiInfoElement element : list) {
            add(element);
        }
    }

    public void add(MeshiInfoElement[] list) {
        for (MeshiInfoElement element : list) {
            add(element);
        }
    }

    public void add(DoubleInfoElementList list) {
        for (DoubleInfoElement element : list) {
            add(element);
        }
    }

    public String toXML(int indentationTabs, MeshiInfoElement header, String endTag) {
        String out = header.toXML() + "\n";
        String tabs = "";
        for (int i = 0; i < indentationTabs; i++)
            tabs += "\t";
        int i = 0;
        for (MeshiInfoElement element : this) {
            if (i % 2 == 0) out += tabs;
            out += (element.toXML() + "                                                                                                       ").substring(0, 80);
            if (i % 2 == 1) out += "\n";
            i++;
        }
        out += endTag;
        return out;
    }

}
