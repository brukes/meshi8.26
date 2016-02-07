/*
 * Copyright (c) 2010. of Chen Keasar, BGU . For free use under LGPL
 */

package meshi.util;

import meshi.energy.EvaluationException;
import meshi.energy.simpleEnergyTerms.angle.AngleParametersList;
import meshi.energy.simpleEnergyTerms.bond.BondParametersList;
import meshi.molecularElements.Protein;
import meshi.sequences.AlignmentException;
import meshi.util.filters.Filter;
import meshi.util.info.ProteinInfo;
import meshi.util.info.ProteinInfoList;

import java.io.File;
import java.io.FileFilter;

/**

 */
public class PdbFilesScanner {
    private final BondParametersList bondParametersList;
    private final File[] files;
    private final AngleParametersList angleParametersList;

    public PdbFilesScanner(File directory, FileFilter filter, CommandList commands) {
        if (!directory.exists()) throw new RuntimeException(directory + " does not a exist.");
        if (!directory.isDirectory()) throw new RuntimeException(directory + " is not a directory.");
        files = directory.listFiles(filter);
        System.out.println("Scanning " + files.length + " files from " + directory.getAbsolutePath());
        if (commands != null) {
            bondParametersList = Utils.getBondParameters(commands);
            angleParametersList = Utils.getAngleParameters(commands);
        } else {
            bondParametersList = null;
            angleParametersList = null;
        }
    }

    public PdbFilesScanner(File directory, FileFilter filter) {
        this(directory, filter, null);
    }

    public ProteinInfoList analyze(ProteinAnalyzer analyzer) throws UpdateableException, EvaluationException,AlignmentException {
        return analyze(analyzer, null);
    }

    public ProteinInfoList analyze(ProteinAnalyzer analyzer, Filter fileFilter) throws UpdateableException,EvaluationException,AlignmentException {
        ProteinInfoList proteinInfoList = analyzer.getProteinInfoList();
        ProteinInfo proteinInfo;
        int i = 0;
        for (File file : files) {
            if ((fileFilter == null) || (fileFilter.accept(file))) {
                //System.out.println("Analyzing "+file.getAbsolutePath());
                System.out.print(".");
                if ((i++) % 50 == 0) System.out.println(" " + i + " ");
                proteinInfo = analyze(file, analyzer);
                if (proteinInfo != null)
                    proteinInfoList.add(proteinInfo);
            } else {
                System.out.println("Ignoring " + file.getAbsolutePath());
            }
        }
        return proteinInfoList;
    }

    public ProteinInfo analyze(File file, ProteinAnalyzer analyzer) throws UpdateableException, EvaluationException,AlignmentException {
        Protein model;
        try {
            model = Protein.getExtendedAtomsProteinFromPdbFile(file, bondParametersList, angleParametersList);
        } catch (Exception ex) {
            System.out.println("Failed to analyze " + file + " due to " + ex);
            ex.printStackTrace();
            return null;
        }
        return analyzer.analyze(model);
    }


}
