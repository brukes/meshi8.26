package meshi.sequences.aligner;

import meshi.util.file.MeshiLineReader;

import java.io.IOException;

public class IdentityMatrix extends SubstitutionMatrix{
	public final static int asciiA=65;
    public static final double GAP_START_PENALTY = -10.0;
    public static final double IN_START_PENALTY = -0.5;


    public IdentityMatrix()  {
        super(GAP_START_PENALTY,IN_START_PENALTY);
        for (int iRow = 0; iRow < substitutionMatrix.length;iRow++)
            for (int jColumn = 0; jColumn < substitutionMatrix.length;jColumn++)
                if (iRow == jColumn) substitutionMatrix[iRow][jColumn]=1;
                else substitutionMatrix[iRow][jColumn]=0;
        String temp = "ARNDCQEGHILKMFPSTWYVBZX";
        order = temp.toCharArray();
        setLetterToIndex(order);
    }
}
