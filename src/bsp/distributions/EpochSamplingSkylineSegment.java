package bsp.distributions;


import beast.math.Binomial;

import javax.swing.text.Segment;

/**
 *
 *
 *
 * Segment is from $(start,end]$ where $end - start = width$
 */
public class EpochSamplingSkylineSegment {

    protected double start,
                     end,
                     width;


    protected int lineageCount;

    protected SegmentType segmentType;

    protected double popSize,
                     samplingIntensity;



    /* Constructors */

    public EpochSamplingSkylineSegment(double end, SegmentType segmentType) {
        this.end = end;
        this.segmentType = segmentType;
        this.start = 0;
        this.width = end;
        this.lineageCount = 0;
    }

    public EpochSamplingSkylineSegment(double start, double end, double width,
                                       int lineageCount, SegmentType segmentType,
                                       double popSize, double samplingIntensity) {
        this.start = start;
        this.end = end;
        this.width = width;
        this.lineageCount = lineageCount;
        this.segmentType = segmentType;
        this.popSize = popSize;
        this.samplingIntensity = samplingIntensity;
    }

    /**
     * Calculate the likelihood of a segment under the linear preferential sampling Kingman coalescent model.
     *
     * @return
     */
    public double calculateLikelihood() {

        final double kchoose2 = Binomial.choose2(lineageCount);

        double lk = -width* (kchoose2/popSize + samplingIntensity*popSize);

        switch (segmentType) {
            case COALESCENT:
                lk += -Math.log(popSize);
                break;

            case SAMPLE:
                lk += Math.log(samplingIntensity*popSize);
                break;

            default:
                break;
        }
        //System.out.println(lk + "-\t" + type+"\t"+lineageCount);

        return lk;
    }


    /* Getters */
    public double getStart() {
        return start;
    }

    public double getEnd() {
        return end;
    }

    public double getWidth() {
        return width;
    }

    public int getLineageCount() {
        return lineageCount;
    }

    public SegmentType getSegmentType() {
        return segmentType;
    }

    public double getPopSize() {
        return popSize;
    }

    public double getSamplingIntensity() {
        return samplingIntensity;
    }


    /* Setters */
    public void setStart(double start) {
        this.start = start;
    }

    public void setEnd(double end) {
        this.end = end;
    }

    public void setWidth(double width) {
        this.width = width;
    }

    public void setLineageCount(int lineageCount) {
        this.lineageCount = lineageCount;
    }

    public void setSegmentType(SegmentType segmentType) {
        this.segmentType = segmentType;
    }

    public void setPopSize(double popSize) {
        this.popSize = popSize;
    }

    public void setSamplingIntensity(double samplingIntensity) {
        this.samplingIntensity = samplingIntensity;
    }


    /* Output */
    public String toString() {
        //return "segment(" + start + "-" + end + ":" + width +" / " + lineageCount + "," + segmentType +")";
        return String.format("%10s%10s%10s |%10d%12s |%10s%18s",
                             start, end, width, lineageCount, segmentType, popSize, samplingIntensity);
    }

}