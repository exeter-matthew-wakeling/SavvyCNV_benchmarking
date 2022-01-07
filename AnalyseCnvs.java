import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

/**
 * Analyse the CNV calls made by SavvyCNV, compared against the mlpa results.
 * Count true/false positive/negatives.
 *
 * @author Matthew Wakeling
 */
public class AnalyseCnvs
{
	public static final Pattern TAB = Pattern.compile("\t");

	public static void main(String[] args) throws Exception {
		BufferedReader in = new BufferedReader("-".equals(args[0]) ? new InputStreamReader(System.in) : new FileReader(args[0]));
		boolean decon = false;
		boolean cnvkit = false;
		boolean excavator2 = false;
		boolean copywriter = false;
		boolean joint = false;
		String truthFile = "truth_set.csv";
		String commonFile = "common_cnvs";
		List<Integer> optimiseColumns = new ArrayList<Integer>();
		for (int i = 1; i < args.length; i++) {
			if ("-decon".equals(args[i])) {
				decon = true;
			} else if ("-cnvkit".equals(args[i])) {
				cnvkit = true;
			} else if ("-excavator2".equals(args[i])) {
				excavator2 = true;
			} else if ("-copywriter".equals(args[i])) {
				copywriter = true;
			} else if ("-joint".equals(args[i])) {
				joint = true;
			} else if ("-truth".equals(args[i])) {
				i++;
				truthFile = args[i];
			} else if ("-optimise".equals(args[i])) {
				i++;
				optimiseColumns.add(Integer.parseInt(args[i]) - 1);
			}
		}
		Map<String, List<Call>> calls = new HashMap<String, List<Call>>();
		String line = in.readLine();
		while (line != null) {
			String[] split = TAB.split(line);
			String sampleName = null;
			Call call = null;
			if (decon) {
				sampleName = split[1];
				if (!"Sample".equals(sampleName)) {
					String chromosome = split[10];
					int start = Integer.parseInt(split[8]);
					int end = Integer.parseInt(split[9]);
					String dupDel = split[6];
					if ("deletion".equals(dupDel)) {
						dupDel = "Deletion";
					} else if ("duplication".equals(dupDel)) {
						dupDel = "Duplication";
					} else {
						throw new Exception("Unrecognised CNV type \"" + dupDel + "\" on line \"" + line + "\"");
					}
					call = new Call(chromosome, start, end, dupDel, line);
				}
			} else if (cnvkit || excavator2 || copywriter) {
				sampleName = split[0];
				String chromosome = split[1];
				int start = Integer.parseInt(split[2]);
				int end = Integer.parseInt(split[3]);
				String dupDel = copywriter ? split[4] : (Double.parseDouble(split[(excavator2 ? 4 : 5)]) > 0.0 ? "Duplication" : "Deletion");
				call = new Call(chromosome, start, end, dupDel, line);
			} else {
				if (joint) {
					sampleName = split[9];
					sampleName = sampleName.substring(sampleName.indexOf("/") + 1);
				} else {
					sampleName = split[9];
				}
				String chromosome = split[0];
				if ("Xmale".equals(chromosome) || "Xfemale".equals(chromosome)) {
					chromosome = "X";
				}
				int start = Integer.parseInt(split[1]);
				int end = Integer.parseInt(split[2]);
				if ("Deletion".equals(split[3]) || "Duplication".equals(split[3])) {
					call = new Call(chromosome, start, end, split[3], line);
				}
			}
			if (!((call == null) || call.getChromosome().startsWith("Y") || call.getChromosome().startsWith("G") || call.getChromosome().startsWith("M"))) {
				List<Call> sampleCalls = calls.get(sampleName);
				if (sampleCalls == null) {
					sampleCalls = new ArrayList<Call>();
					calls.put(sampleName, sampleCalls);
				}
				if (optimiseColumns != null) {
					double[] columns = new double[optimiseColumns.size()];
					for (int i = 0; i < columns.length; i++) {
						if (optimiseColumns.get(i) == -1) {
							if (cnvkit) {
								columns[i] = Double.parseDouble(split[9]) / (Double.parseDouble(split[3]) - Double.parseDouble(split[2]));
							} else if (excavator2) {
								columns[i] = Double.parseDouble(split[3]) - Double.parseDouble(split[2]);
							} else if (joint) {
								columns[i] = Math.max(Double.parseDouble(split[5]), Double.parseDouble(split[6]));
							} else {
								columns[i] = Double.parseDouble(split[12]) / (Double.parseDouble(split[9]) - Double.parseDouble(split[8]));
							}
						} else {
							columns[i] = Double.parseDouble(split[optimiseColumns.get(i)]);
						}
					}
					call.setOptimiseColumns(columns);
				}
				sampleCalls.add(call);
			}
			line = in.readLine();
		}
		Map<String, List<Call>> trueCalls = new HashMap<String, List<Call>>();
		in = new BufferedReader(new FileReader(truthFile));
		line = in.readLine();
		while (line != null) {
			String[] split = TAB.split(line);
			if ("Normal".equals(split[4]) || "ExonCNV".equals(split[4])) {
				String sampleName = split[0];
				String chr = split[7];
				int start = Integer.parseInt(split[8]);
				int end = Integer.parseInt(split[9]);
				String dupDel = split[5];
				List<Call> trueSampleCalls = trueCalls.get(sampleName);
				if (trueSampleCalls == null) {
					trueSampleCalls = new ArrayList<Call>();
					trueCalls.put(sampleName, trueSampleCalls);
				}
				trueSampleCalls.add(new Call(chr, start, end, dupDel, line));
			}
			line = in.readLine();
		}
		in.close();
		in = new BufferedReader(new FileReader(commonFile));
		line = in.readLine();
		List<Call> commonCalls = new ArrayList<Call>();
		while (line != null) {
			String[] split = TAB.split(line);
			commonCalls.add(new Call(split[0], Integer.parseInt(split[1]), Integer.parseInt(split[2]), "", line));
			line = in.readLine();
		}
		in.close();
		List<Call> fp = new ArrayList<Call>();
		Map<Call, Set<Call>> tp = new HashMap<Call, Set<Call>>();
		for (String sampleName : calls.keySet()) {
			List<Call> sampleCalls = calls.get(sampleName);
			List<Call> trueSampleCalls = trueCalls.get(sampleName);
			if (trueSampleCalls == null) {
				trueSampleCalls = new ArrayList<Call>();
			}
			for (Call call : sampleCalls) {
				boolean common = false;
				for (Call commonCall : commonCalls) {
					if (call.containedIn(commonCall)) {
						common = true;
					}
				}
				boolean overlaps = false;
				Set<Call> truths = new HashSet<Call>();
				for (Call trueCall : trueSampleCalls) {
					if (call.overlaps(trueCall)) {
						overlaps = true;
						truths.add(trueCall);
					}
				}
				if (overlaps) {
					tp.put(call, truths);
				} else if (!common) {
					fp.add(call);
				}
			}
		}
		Set<Call> truths = new HashSet<Call>();
		for (Set<Call> callTruths : tp.values()) {
			truths.addAll(callTruths);
		}
		if (optimiseColumns.isEmpty()) {
			System.out.println(truths.size() + "\t" + fp.size());
		} else {
			int[] bestFp = new int[truths.size() + 1];
			for (int i = 0; i < bestFp.length; i++) {
				bestFp[i] = Integer.MAX_VALUE;
			}
			double[][] bestParams = new double[truths.size() + 1][optimiseColumns.size()];
			optimise(tp, fp, new ArrayList<Double>(), bestFp, bestParams);
			for (int i = 1; i < bestFp.length; i++) {
				System.out.print(i + "\t" + bestFp[i]);
				for (int o = 0; o < bestParams[i].length; o++) {
					System.out.print("\t" + bestParams[i][o]);
				}
				System.out.println("");
			}
		}
	}

	public static class Call
	{
		String chr, dupDel, line;
		int start, end;
		double[] optimiseColumns;

		public Call(String chr, int start, int end, String dupDel, String line) {
			this.chr = chr;
			this.start = start;
			this.end = end;
			this.dupDel = dupDel;
			this.line = line;
		}

		public String getChromosome() {
			return chr;
		}

		public int getLength() {
			return end - start;
		}

		public String getLine() {
			return line;
		}

		public void setOptimiseColumns(double[] optimiseColumns) {
			this.optimiseColumns = optimiseColumns;
		}

		public double[] getOptimiseColumns() {
			return optimiseColumns;
		}

		public boolean overlaps(String oChr, int oStart, int oEnd, String oDupDel) {
			return chr.equals(oChr) && (oStart <= end) && (start <= oEnd) && ("".equals(oDupDel) || "".equals(dupDel) || dupDel.toLowerCase().equals(oDupDel.toLowerCase()));
		}

		public boolean overlaps(Call o) {
			return overlaps(o.chr, o.start, o.end, o.dupDel);
		}

		public boolean containedIn(Call o) {
			return chr.equals(o.chr) && (o.start <= start) && (end <= o.end) && ("".equals(o.dupDel) || "".equals(dupDel) || dupDel.toLowerCase().equals(o.dupDel.toLowerCase()));
		}

		public String toString() {
			return chr + ":" + start + "-" + end + "(" + dupDel + ")";
		}
	}

	public static void optimise(Map<Call, Set<Call>> tp, List<Call> fp, List<Double> optimised, int[] bestFp, double[][] bestParams) {
		int nextColumn = optimised.size();
		TreeSet<Double> values = new TreeSet<Double>();
		for (Call call : tp.keySet()) {
			values.add(call.getOptimiseColumns()[nextColumn]);
		}
		for (Double threshold : values) {
			Map<Call, Set<Call>> newTp = new HashMap<Call, Set<Call>>();
			for (Map.Entry<Call, Set<Call>> entry : tp.entrySet()) {
				if (entry.getKey().getOptimiseColumns()[nextColumn] >= threshold) {
					newTp.put(entry.getKey(), entry.getValue());
				}
			}
			List<Call> newFp = new ArrayList<Call>();
			for (Call call : fp) {
				if (call.getOptimiseColumns()[nextColumn] >= threshold) {
					newFp.add(call);
				}
			}
			List<Double> newOptimised = new ArrayList<Double>();
			newOptimised.addAll(optimised);
			newOptimised.add(threshold);
			if (optimised.size() + 1 < bestParams[0].length) {
				optimise(newTp, newFp, newOptimised, bestFp, bestParams);
			} else {
				Set<Call> truths = new HashSet<Call>();
				for (Set<Call> callTruths : newTp.values()) {
					truths.addAll(callTruths);
				}
				for (int t = 0; t <= truths.size(); t++) {
					if (bestFp[t] > newFp.size()) {
						bestFp[t] = newFp.size();
						for (int i = 0; i < newOptimised.size(); i++) {
							bestParams[t][i] = newOptimised.get(i);
						}
					}
				}
			}
		}
	}
}
