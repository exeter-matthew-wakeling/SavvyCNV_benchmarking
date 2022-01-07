import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.List;
import java.util.regex.Pattern;

/**
 * Process the output from GATK GermlineCNVCaller to produce CNV calls.
 *
 * @author Matthew Wakeling
 */
public class ProcessGatkGcnvResults
{
	public static final Pattern TAB = Pattern.compile("\t");

	public static void main(String[] args) throws Exception {
		// First argument is a name of a directory containing results.
		// Second argument is male/female
		boolean male = "male".equals(args[1]);
		File baseDir = new File(args[0]);
		List<Interval<String>> intervals = new ArrayList<Interval<String>>();
		BufferedReader in = new BufferedReader(new FileReader(args[0] + "/interval_list.tsv"));
		String line = in.readLine();
		while (line != null) {
			if (!(line.startsWith("@") || line.startsWith("CONTIG"))) {
				String[] split = TAB.split(line);
				intervals.add(new Interval<String>(split[0], Integer.parseInt(split[1]), Integer.parseInt(split[2]), null));
			}
			line = in.readLine();
		}
		in.close();
		File[] dirList = baseDir.listFiles();
		for (File dirEntry : dirList) {
			if (dirEntry.isDirectory() && dirEntry.getName().startsWith("SAMPLE_")) {
				in = new BufferedReader(new FileReader(dirEntry.getPath() + "/sample_name.txt"));
				String sampleName = in.readLine();
				System.err.println("Found sample \"" + sampleName + "\"");
				in.close();
				in = new BufferedReader(new FileReader(dirEntry.getPath() + "/baseline_copy_number_t.tsv"));
				in.readLine();
				in.readLine();
				int[] baseline = new int[intervals.size()];
				for (int i = 0; i < baseline.length; i++) {
					baseline[i] = Integer.parseInt(in.readLine());
				}
				in.close();
				in = new BufferedReader(new FileReader(dirEntry.getPath() + "/log_c_emission_tc.tsv"));
				in.readLine();
				in.readLine();
				int state = 0;
				// -1 for deletion, 0 for normal, 1 for insertion.
				String lastChr = "";
				int startCnv = -1;
				int endCnv = -1;
				double score = 0.0;
				int buckets = 0;
				for (int i = 0; i < intervals.size(); i++) {
					line = in.readLine();
					String[] split = TAB.split(line);
					double[] gqs = new double[6];
					for (int o = 0; o < 6; o++) {
						gqs[o] = Double.parseDouble(split[o]);
					}
					int normalCN = 2;
					String chr = intervals.get(i).getChromosome();
					if (!chr.equals(lastChr)) {
						if (state < 0) {
							System.out.println(lastChr + "\t" + startCnv + "\t" + endCnv + "\tDeletion\t" + score + "\t" + buckets + "\t" + (endCnv - startCnv) + "\t" + (score / buckets) + "\t" + (score * 1000000.0 / (endCnv - startCnv)) + "\t" + sampleName);
						} else if (state > 0) {
							System.out.println(lastChr + "\t" + startCnv + "\t" + endCnv + "\tDuplication\t" + score + "\t" + buckets + "\t" + (endCnv - startCnv) + "\t" + (score / buckets) + "\t" + (score * 1000000.0 / (endCnv - startCnv)) + "\t" + sampleName);
						}
						lastChr = chr;
						state = 0;
						score = 0.0;
						buckets = 0;
						startCnv = -1;
						endCnv = -1;
					}
					if (male) {
						if ("X".equals(chr) || "Y".equals(chr)) {
							normalCN = 1;
						}
					}
					if ((!"Y".equals(chr)) || male) {
						double delScore = -10000.0;
						double normScore = -10000.0;
						double dupScore = -10000.0;
						for (int o = 0; o < 6; o++) {
							if (o < normalCN) {
								delScore = Math.max(delScore, gqs[o]);
							} else if (o > normalCN) {
								dupScore = Math.max(dupScore, gqs[o]);
							} else {
								normScore = Math.max(normScore, gqs[o]);
							}
						}
						if (state < 0) {
							if ((delScore > normScore) && (delScore > dupScore)) {
								// Still a deletion
								score += (delScore - Math.max(normScore, dupScore));
								buckets++;
								endCnv = intervals.get(i).getEnd();
							} else if ((normScore > delScore) && (normScore > dupScore)) {
								// Transition to normal. End CNV
								System.out.println(chr + "\t" + startCnv + "\t" + endCnv + "\tDeletion\t" + score + "\t" + buckets + "\t" + (endCnv - startCnv) + "\t" + (score / buckets) + "\t" + (score * 1000000.0 / (endCnv - startCnv)) + "\t" + sampleName);
								state = 0;
								score = 0.0;
								buckets = 0;
								startCnv = -1;
								endCnv = -1;
							} else {
								// Transition to duplication. End CNV
								System.out.println(chr + "\t" + startCnv + "\t" + endCnv + "\tDeletion\t" + score + "\t" + buckets + "\t" + (endCnv - startCnv) + "\t" + (score / buckets) + "\t" + (score * 1000000.0 / (endCnv - startCnv)) + "\t" + sampleName);
								state = 1;
								score = (dupScore - Math.max(normScore, delScore));
								buckets = 1;
								startCnv = intervals.get(i).getStart() - 1;
								endCnv = intervals.get(i).getEnd();
							}
						} else if (state > 0) {
							if ((delScore > normScore) && (delScore > dupScore)) {
								// Transition to deletion
								System.out.println(chr + "\t" + startCnv + "\t" + endCnv + "\tDuplication\t" + score + "\t" + buckets + "\t" + (endCnv - startCnv) + "\t" + (score / buckets) + "\t" + (score * 1000000.0 / (endCnv - startCnv)) + "\t" + sampleName);
								state = -1;
								score = (delScore - Math.max(normScore, dupScore));
								buckets = 1;
								startCnv = intervals.get(i).getStart() - 1;
								endCnv = intervals.get(i).getEnd();
							} else if ((normScore > delScore) && (normScore > dupScore)) {
								// Transition to normal. End CNV
								System.out.println(chr + "\t" + startCnv + "\t" + endCnv + "\tDuplication\t" + score + "\t" + buckets + "\t" + (endCnv - startCnv) + "\t" + (score / buckets) + "\t" + (score * 1000000.0 / (endCnv - startCnv)) + "\t" + sampleName);
								state = 0;
								score = 0.0;
								buckets = 0;
								startCnv = -1;
								endCnv = -1;
							} else {
								// Still a duplication
								score += (dupScore - Math.max(normScore, delScore));
								buckets++;
								endCnv = intervals.get(i).getEnd();
							}
						} else {
							if ((delScore > normScore) && (delScore > dupScore)) {
								// Transition to deletion
								state = -1;
								score = (delScore - Math.max(normScore, dupScore));
								buckets = 1;
								startCnv = intervals.get(i).getStart() - 1;
								endCnv = intervals.get(i).getEnd();
							} else if ((dupScore > delScore) && (dupScore > normScore)) {
								state = 1;
								score = (dupScore - Math.max(normScore, delScore));
								buckets = 1;
								startCnv = intervals.get(i).getStart() - 1;
								endCnv = intervals.get(i).getEnd();
							}
						}
					}
				}
				in.close();
				if (state < 0) {
					System.out.println(lastChr + "\t" + startCnv + "\t" + endCnv + "\tDeletion\t" + score + "\t" + buckets + "\t" + (endCnv - startCnv) + "\t" + (score / buckets) + "\t" + (score * 1000000.0 / (endCnv - startCnv)) + "\t" + sampleName);
				} else if (state > 0) {
					System.out.println(lastChr + "\t" + startCnv + "\t" + endCnv + "\tDuplication\t" + score + "\t" + buckets + "\t" + (endCnv - startCnv) + "\t" + (score / buckets) + "\t" + (score * 1000000.0 / (endCnv - startCnv)) + "\t" + sampleName);
				}
			}
		}
	}
}
