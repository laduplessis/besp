package bsp.distributions;

import beast.math.Binomial;

/**
 * Tree Interval (segment) with associated parameterisation for the BSP
 *
 * @author Louis du Plessis
 * @date 2019/01/17
 */
public class BSPSegment {

    protected double start,
            end,
            width;

    protected int lineageCount;

    protected SegmentType segmentType;

    protected double popSize;


    /* Constructors */
    public BSPSegment() {

    }

    public BSPSegment(double end, SegmentType segmentType) {
        this.end = end;
        this.segmentType = segmentType;
        this.start = 0;
        this.width = end;
        this.lineageCount = 0;
    }

    public BSPSegment(double start, double end, double width,
                                       int lineageCount, SegmentType segmentType,
                                       double popSize) {
        this.start = start;
        this.end = end;
        this.width = width;
        this.lineageCount = lineageCount;
        this.segmentType = segmentType;
        this.popSize = popSize;
    }

    public BSPSegment(BSPSegment segment) {
        this.start        = segment.start;
        this.end          = segment.end;
        this.width        = segment.width;
        this.lineageCount = segment.lineageCount;
        this.segmentType  = segment.segmentType;
        this.popSize      = segment.popSize;
    }


    /**
     * Calculate the likelihood of a segment under a fixed population size Kingman coalescent model.
     *
     * @return
     */
    public double calculateLikelihood() {

        final double kchoose2 = Binomial.choose2(lineageCount);

        double lk = -width* (kchoose2/popSize);

        switch (segmentType) {
            case COALESCENT:
                lk += -Math.log(popSize);
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


    /* Output */
    public String toString() {
        //return "segment(" + start + "-" + end + ":" + width +" / " + lineageCount + "," + segmentType +")";
        return String.format("%10s%10s%10s |%10d%12s",
                start, end, width, lineageCount, segmentType, popSize);
    }


}
