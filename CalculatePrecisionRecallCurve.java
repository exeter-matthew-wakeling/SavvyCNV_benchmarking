import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.regex.Pattern;

/**
 * Calculates the best precision-recall curve that can be constructed from a set of separate tests.
 * The input is read from the standard input, and is expected to be tab-separated lines.
 * Each line is a separate test. The number of true positive detections and false positives are read from the line.
 * The first argument is the column number of the true positives.
 * The second argument is the column number of the false positives.
 * The third argument is total number of condition positives (true positives plus false negatives).
 * The fourth argument is optional - if it is not present then the software will emit just the input lines that are on the best precision-recall curve.
 * If the fourth argument is present, then we recommend a value of 1, and the software will interpolate lines between these points.
 * The lines will be curved on a precision-recall graph, but would be straight lines on a ROC curve (but a ROC curve cannot be created without having a well-defined number of true negatives).
 *
 * @author Matthew Wakeling
 */
public class CalculatePrecisionRecallCurve
{
	public static final Pattern TAB = Pattern.compile("\t");

	public static void main(String[] args) throws Exception {
		int truePosColumn = Integer.parseInt(args[0]) - 1;
		int falsePosColumn = Integer.parseInt(args[1]) - 1;
		int totalPos = Integer.parseInt(args[2]);
		int intermediate = 0;
		if (args.length > 3) {
			intermediate = Integer.parseInt(args[3]);
		}
		ArrayList<Point> points = new ArrayList<Point>();
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		String line = in.readLine();
		while (line != null) {
			String[] split = TAB.split(line);
			if (split.length > Math.max(truePosColumn, falsePosColumn)) {
				points.add(new Point(Integer.parseInt(split[truePosColumn]), Integer.parseInt(split[falsePosColumn]), totalPos, line));
			}
			line = in.readLine();
		}
		double currentX = 0.0;
		double currentY = 0.0;
		Point lastPoint = new Point(0, 0, totalPos, "");
		while (currentY < 1.0) {
			double bestX = 1.0;
			double bestY = 1.0;
			Point bestPoint = new Point(totalPos, Integer.MAX_VALUE, totalPos, "");
			boolean pointFound = false;
			for (Point p : points) {
				if ((p.getX() >= currentX) && (p.getY() > currentY) && ((p.getY() > bestY) || (bestX > currentX)) && ((p.getY() - currentY) / (p.getX() - currentX) >= (bestY - currentY) / (bestX - currentX))) {
					bestX = p.getX();
					bestY = p.getY();
					bestPoint = p;
					pointFound = true;
				}
			}
			if (pointFound) {
				if (intermediate > 0) {
					int intermediate2 = intermediate * Math.min(100, Math.max(bestPoint.getTruePos() - lastPoint.getTruePos(), bestPoint.getFalsePos() - lastPoint.getFalsePos()));
					for (int i = 1; i < intermediate2; i++) {
						double fraction = (1.0 * i) / intermediate2;
						double truePos = lastPoint.getTruePos() + fraction * (bestPoint.getTruePos() - lastPoint.getTruePos());
						double falseNeg = lastPoint.getFalseNeg() + fraction * (bestPoint.getFalseNeg() - lastPoint.getFalseNeg());
						double falsePos = lastPoint.getFalsePos() + fraction * (bestPoint.getFalsePos() - lastPoint.getFalsePos());
						System.out.println((truePos / (truePos + falsePos)) + "\t" + (truePos / (truePos + falseNeg)));
					}
				}
				System.out.println(bestPoint.getPrecision() + "\t" + bestPoint.getY() + "\t" + bestPoint.getLine());
			}
			currentX = bestX;
			currentY = bestY;
			lastPoint = bestPoint;
		}
	}

	public static class Point
	{
		private int truePos, falsePos, falseNeg, trueNeg;
		private double x, y;
		private String line;

		public Point(int truePos, int falsePos, int totalPos, String line) {
			this.truePos = truePos;
			this.falsePos = falsePos;
			this.falseNeg = totalPos - truePos;
			this.trueNeg = 10000000 - falsePos;
			this.line = line;
		}

		public double getX() {
			return (1.0 * falsePos) / (falsePos + trueNeg);
		}

		public double getY() {
			return (1.0 * truePos) / (truePos + falseNeg);
		}

		public double getPrecision() {
			return (1.0 * truePos) / (truePos + falsePos);
		}

		public int getTruePos() {
			return truePos;
		}

		public int getFalsePos() {
			return falsePos;
		}

		public int getFalseNeg() {
			return falseNeg;
		}

		public String getLine() {
			return line;
		}

		public String toString() {
			return "(" + getX() + ", " + getY() + ")";
		}
	}
}
