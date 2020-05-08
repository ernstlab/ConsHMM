import edu.mit.compbio.ChromHMM.*;
import java.util.*;
import java.io.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class ReassignVariantState
{

//static int[][] data = {
//	{0,0,0,0,0,0},
//	{0,0,1,0,0,0},
//	{0,0,1,0,0,0},
//	{0,0,1,0,0,0},
//	{0,0,0,0,0,0},
//	{0,0,0,0,0,0},
//{0,0,0,0,0,0}};

    static String[] NUCLEOTIDES = {"A", "C", "T", "G"};

    static LinkedList<int[]> currentdata = new LinkedList<int[]> ();
    static LinkedList<String[]> currentseq = new LinkedList<String[]> ();
    static int posReference; // will keep track of position in line of reference species

    static int[] createFeatures(String[] splitLine, int posReference, String nucleotide) {
        int[] featureArray = new int[(splitLine.length - 3) * 2];
        int pos = 0;
        
        for (int i = 2; i < splitLine.length; i++) {
            if (i != posReference) {
                if ((splitLine[i].compareTo("-") == 0) || (splitLine[i].compareTo("X") == 0)) { // species does not align at all
                    featureArray[pos] = 0;
                    featureArray[pos + 1] = 2;
                } else {
                    featureArray[pos] = 1;
                    if (splitLine[i].compareTo(nucleotide) == 0) {
                        featureArray[pos + 1] = 1;
                    } else {
                        featureArray[pos + 1] = 0;
                    }
                }
                pos += 2;
            }
        }

        return featureArray;
    }

    static void writeToFile(GZIPOutputStream out, int pos, ChromHMM theChromHMM, int windowsize) throws IOException {
        out.write(((pos - (windowsize / 2)) + " ").getBytes()); 
        for (int j = 0; j < 4; j++) {
            String nucleotide = NUCLEOTIDES[j];
            int[] newDataLine = createFeatures(currentseq.get(windowsize / 2), posReference, nucleotide);
            currentdata.set(windowsize / 2, newDataLine);
        
            int[][] arrayData = currentdata.toArray(new int[currentdata.size()][]);
            // System.out.println(arrayData[0][0]); 
            int nmaxstate = theChromHMM.getMaxStateAtPos(arrayData, windowsize / 2);
    
            out.write((nmaxstate + " ").getBytes());
        }
        out.write("\n".getBytes());
    }

    static int[] getChromSizes(String sizesfile, String chrom, int chunksize) throws IOException {
        int[] chromInfo = new int[2]; // chromInfo[0] = chromosome size; chromInfo[1] = total number of chunks
        try (BufferedReader bufferedReader = new BufferedReader(new FileReader(sizesfile))) {
            String line = bufferedReader.readLine();
            while (line != null) {
                String[] splitLine = line.trim().split("\t", -1);
                String curchrom = splitLine[0];
                int chromsize = Integer.parseInt(splitLine[1]); // chromosome size
                int numChunks = chromsize / chunksize + 1; // total number of chunks for this chromosome

                if (curchrom.compareTo(chrom) == 0) {
                    chromInfo[0] = chromsize;
                    chromInfo[1] = numChunks;
                    return chromInfo;
                }
                line = bufferedReader.readLine();
            }
        } catch (FileNotFoundException e) {
            System.out.println(sizesfile + " not found.");
        } catch (IOException e) {
            System.out.println("Couldn't read from file.");
        }
        chromInfo[0] = -1;
        chromInfo[1] = -1;
        return chromInfo;
    }

    public static void main(String[] args) throws IOException {
        if (args.length != 9) {
            System.out.println("Usage: reassignVariants <input file> <output file> <reference species> <model file name> <chunk index> <chunk size> <chrom size file> <current chromosome> <window size>");
            System.exit(0);
        }
        
        String inputfile = args[0];
        String outputfile = args[1];
        String reference = args[2];
        String modelfile = args[3];
        int chunk = Integer.parseInt(args[4]); // needed because we are splitting up the MAF sequence files
        int chunksize = Integer.parseInt(args[5]);
        String sizesfile = args[6];
        String curchrom = args[7];
        int windowsize = Integer.parseInt(args[8]);

        int[] chromInfo = getChromSizes(sizesfile, curchrom, chunksize);

        // Load ChromHMM model	
	    ChromHMM theChromHMM = new ChromHMM(modelfile);
           
        // Set up reading input file
        GZIPInputStream in = new GZIPInputStream(new FileInputStream(inputfile));
        GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(outputfile));
        out.write("pos A C T G\n".getBytes());

        Reader decoder = new InputStreamReader(in);
        BufferedReader br = new BufferedReader(decoder);

        String line;        
        int lastpos = -11;
        String[] lastLine;
        int last_length = 0;
        String last_chrom = "chr";

        while ((line = br.readLine()) != null) {
            String[] splitLine = line.trim().split(",", -1);
            lastLine = splitLine;
            last_length = lastLine.length;
            last_chrom = lastLine[0];
            // Test if it's a header line
            try {
                int test = Integer.parseInt(splitLine[1]);
            } catch (NumberFormatException | NullPointerException nfe) {
                for (int j = 0; j < splitLine.length; j++) {
                    if (splitLine[j].compareTo(reference) == 0) {
                        posReference = j;
                    }
                }
                continue;
            }
            int curpos = Integer.parseInt(splitLine[1]);

            // If there is a gap fill in those positions with empty lines
            
            // First need to test if this is the first chromosome or not
            if (chunk != 0) {
                if (lastpos < 0) { // don't need to fill gaps before the first position if this isn't the first chunk
                    lastpos = curpos - 11; // setting this so the next if statement will always be false
                }
            }

            if (curpos != lastpos + 1) {
                String[] emptyLine = new String[splitLine.length];
                for (int i = lastpos + 1; i < curpos; i++) {
                    emptyLine[0] = splitLine[0]; // copy chromosome
                    emptyLine[1] = Integer.toString(i); // fill in position
                    for (int j = 2; j < emptyLine.length; j++) {
                        emptyLine[j] = "-";
                    }

                    int[] dataLine = createFeatures(emptyLine, posReference, "N");
                    currentdata.addLast(dataLine);
                    currentseq.addLast(emptyLine);

                    if (currentdata.size() < windowsize) {
                        continue;
                    }

                    writeToFile(out, i, theChromHMM, windowsize);
                    
                    currentdata.pop();
                    currentseq.pop();
                }
            }
            //System.out.println("b:");
            
            // If done with empty lines process current line 
            int[] dataLine = createFeatures(splitLine, posReference, splitLine[posReference]);

            //System.out.println("c:");
            currentdata.addLast(dataLine); // add to end of LL
            currentseq.addLast(splitLine);
            lastpos = curpos;

            if (currentdata.size() < windowsize) { // keep adding lines if we're below window size
                continue;
            }

            writeToFile(out, curpos, theChromHMM, windowsize);
            //System.out.println("d:");
            // Assign state to position in the middle of window
            
            currentdata.pop(); // remove first element because a new one will be added during the next iteration
            currentseq.pop();
            //System.out.println("dd");
        }
        
        //System.out.println("e:");
        // Add additional WINDOWSIZE / 2 empty lines at the end so that the last positions get printed
        String[] emptyLine = new String[last_length];
        for (int i = lastpos + 1; i < lastpos + (windowsize / 2) + 1; i++) {
            emptyLine[0] = last_chrom; // copy chromosome
            emptyLine[1] = Integer.toString(i); // fill in position
            for (int j = 2; j < emptyLine.length; j++) {
                emptyLine[j] = "-";
            }

            int[] dataLine = createFeatures(emptyLine, posReference, "N");
            currentdata.addLast(dataLine);
            currentseq.addLast(emptyLine);
            
            writeToFile(out, i, theChromHMM, windowsize);
            
            currentdata.pop();
            currentseq.pop();
        }
        //System.out.println("f:");
        // If this is the last chunk, fill with empty lines until the end of the chromosome
        if (chunk == chromInfo[1] - 1) {
            if (lastpos != chromInfo[0] - 1) { // if the last position printed isn't the end of the chromosome
                for (int i = lastpos + (windowsize / 2); i < chromInfo[0] + (windowsize / 2); i++) { // keep addding empty lines
                    emptyLine[0] = last_chrom; // copy chromosome
                    emptyLine[1] = Integer.toString(i); // fill in position
                    for (int j = 2; j < emptyLine.length; j++) {
                        emptyLine[j] = "-";
                    }

                    int[] dataLine = createFeatures(emptyLine, posReference, "N");
                    currentdata.addLast(dataLine);
                    currentseq.addLast(emptyLine);
                    
                    writeToFile(out, i, theChromHMM, windowsize);
                    
                    currentdata.pop();
                    currentseq.pop();
                }
            }
        }
        //System.out.println("g:");

        out.finish();
    }
}
