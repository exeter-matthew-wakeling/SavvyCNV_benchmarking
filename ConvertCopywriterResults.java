import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashSet;
import java.util.regex.Pattern;

/**
 * Convert the output data from CopywriteR into a format that can be analysed.
 *
 * @author Matthew Wakeling
 */
public class ConvertCopywriterResults
{
	public static final Pattern SPACES = Pattern.compile(" +");

	public static void main(String[] args) throws Exception {
		HashSet<String> maleSamples = new HashSet<String>();
		{
			BufferedReader in = new BufferedReader(new FileReader("male_samples"));
			String line = in.readLine();
			while (line != null) {
				maleSamples.add(line);
				line = in.readLine();
			}
			in.close();
		}
		int binSize = Integer.parseInt(args[0]);
		File analysisDir = new File("analysis");
		String[] samples = analysisDir.list();
		for (String sampleName : samples) {
			BufferedReader in = new BufferedReader(new FileReader("analysis/" + sampleName + "/results_" + binSize));
			String line = in.readLine();
			while (line != null) {
				String[] split = SPACES.split(line);
				if ((split.length == 7) && (split[1].endsWith(".vs.none"))) {
					String chr = split[2];
					int start = (int) Double.parseDouble(split[3]);
					int end = (int) Double.parseDouble(split[4]);
					int buckets = Integer.parseInt(split[5]);
					double logDosage = Double.parseDouble(split[6]);
					double origDosage = logDosage;
					if ("23".equals(chr)) {
						chr = "X";
						if (maleSamples.contains(sampleName)) {
							logDosage += 1.0;
						}
					} else if ("24".equals(chr)) {
						chr = "Y";
					}
//					if (!("Y".equals(chr))) {
						System.out.println(sampleName + "\t" + chr + "\t" + start + "\t" + end + "\t" + (logDosage < 0 ? "Deletion" : "Duplication") + "\t" + buckets + "\t" + (((end - start) / binSize + 200) / 1000) + "\t" + origDosage + "\t" + Math.abs(logDosage) + "\t" + (Math.abs(logDosage) * buckets));
//					}
				}
				line = in.readLine();
			}
		}
	}
}
