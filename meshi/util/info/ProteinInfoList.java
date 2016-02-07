/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util.info;

import meshi.util.file.MeshiWriter;

import java.io.IOException;
import java.util.ArrayList;

/**

 */
public class ProteinInfoList extends ArrayList<ProteinInfo> {
    public final String name;

    public ProteinInfoList(String name) {
        this.name = name;
    }

    public void print() throws IOException {
        MeshiWriter writer = new MeshiWriter(name + ".pil");
        print(writer);
        writer.close();
    }

    public void filterByOriginalModel() {
        ProteinInfo infoI, infoJ;
        for (int i = 0; i < size(); i++) {
            infoI = get(i);
            if (infoI.on()) {
                for (int j = i + 1; j < size(); j++) {
                    infoJ = get(j);
                    if (infoJ.on()) {
                        if (infoI.originalModel().equals(infoJ.originalModel()))
                            infoJ.turnOff();
                    }
                }
            }
        }
    }

    public void sort() {
        ProteinInfo temp;
        double sumI, sumJ;
        ProteinInfo[] sorted = new ProteinInfo[size()];
        ProteinInfo infoI, infoJ;
        for (int i = 0; i < size(); i++) {
            sorted[i] = get(i);
        }

        for (int i = 0; i < size(); i++) {
            for (int j = 0; j < i; j++) {
                if (sorted[i].sum() > sorted[j].sum()) {
                    temp = sorted[i];
                    sorted[i] = sorted[j];
                    sorted[j] = temp;
                }
            }
        }
        clear();
        for (ProteinInfo info : sorted)
            add(info);
    }

    public void print(MeshiWriter writer) throws IOException {
        if (size() == 0) {
            System.out.println("Empty list " + this);
            return;
        }
        String header = get(0).header();
        writer.println(header);
        for (ProteinInfo proteinInfo : this) {
            if (proteinInfo.on()) {
                if (!proteinInfo.header().equals(header)) {
                    header = proteinInfo.header();
                    writer.println("%" + header);
                }
                writer.println(proteinInfo.values());
                writer.flush();
            }
        }
    }

    public void print(MeshiInfoXMLwriter writer) throws IOException {
        print(writer, 0);
    }

    public void print(MeshiInfoXMLwriter writer, int indentationTabs) throws IOException {
        if (size() == 0) {
            System.out.println("Empty list " + this);
            return;
        }
        writer.println("<ProteinInfoList name=\"" + name + "\">");
        for (ProteinInfo proteinInfo : this) {
            if (proteinInfo.on()) {
                writer.println(proteinInfo.toXML(indentationTabs + 1));
            }
        }
        writer.println("</ProteinInfoList>");
        writer.flush();
    }

    public ProteinInfoList toZscore() {
        ProteinInfo sum = new ProteinInfo(this.get(0), "sum", null, "sum of protein info");
        sum.toZero();
        ProteinInfo sum2 = new ProteinInfo(this.get(0), "sum2", null, "sum of protein info squared");
        sum2.toZero();
        for (ProteinInfo proteinInfo : this) {
            sum.addElementValue(proteinInfo);
            sum2.addElementValue2(proteinInfo);
        }

        ProteinInfo avg = new ProteinInfo(sum, "avg", sum.fileName(), "averaged protein info");
        avg.divideBy(size());

        ProteinInfo avg2 = new ProteinInfo(sum2, "avg2", sum2.fileName(), "averaged protein info squared");
        avg2.divideBy(size());

        ProteinInfo avgavg = new ProteinInfo(avg, "avgavg", avg.fileName(), "squared average of protein info");
        avgavg.square();

        ProteinInfo variance = new ProteinInfo(avg2, "variance", avg2.fileName(), "Variance of protein info.");
        variance.subtract(avgavg);

        ProteinInfo std = new ProteinInfo(variance, "std", variance.fileName(), "Standard deviation of protein info.");
        std.sqrt();

        ProteinInfoList zScores = new ProteinInfoList(name + ".zScores");
        for (ProteinInfo proteinInfo : this) {
            ProteinInfo zScore = new ProteinInfo(proteinInfo, "temp1", proteinInfo.fileName(), "temp2");
            zScore.subtract(avg);
            zScore.divideBy(std);
            zScores.add(zScore);
        }
        return zScores;
    }
}
