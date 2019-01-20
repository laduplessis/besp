package bsp.util;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.HeapSort;
import bsp.distributions.SegmentType;

public class LightweightTreeIntervals extends CalculationNode {

    final public Input<Tree> treeInput = new Input<>("tree", "tree for which to calculate the intervals", Input.Validate.REQUIRED);
    
    protected double      [] intervalTimes;
    protected double      [] intervalWidths;
    protected int         [] lineageCounts;
    protected SegmentType [] intervalTypes;

    protected double      [] storedIntervalTimes;
    protected double      [] storedIntervalWidths;
    protected int         [] storedLineageCounts;
    protected SegmentType [] storedIntervalTypes;

    private boolean intervalsKnown = false;

    public LightweightTreeIntervals() {
        super();
    }

    public LightweightTreeIntervals(Tree tree) {
        init(tree);
    }
    
    @Override
    public void initAndValidate() {
        Tree tree = treeInput.get();
        int nrEvents = tree.getNodeCount();

        intervalTimes  = new double[nrEvents];
        intervalWidths = new double[nrEvents];
        lineageCounts  = new int[nrEvents];
        intervalTypes  = new SegmentType[nrEvents];

        storedIntervalTimes  = new double[nrEvents];
        storedIntervalWidths = new double[nrEvents];
        storedLineageCounts  = new int[nrEvents];
        storedIntervalTypes  = new SegmentType[nrEvents];
        
        updateIntervals();
    }

    public int getIntervalCount() {
        return intervalTimes.length;
    }

    public int getCoalCount() {
        return treeInput.get().getInternalNodeCount();
    }

    public int getLeafCount() {
        return treeInput.get().getLeafNodeCount();
    }

    public double getIntervalWidth(int i) {
        return intervalWidths[i];
    }

    public int getLineageCount(int i) {
        return lineageCounts[i];
    }

    public SegmentType getIntervalType(int i) {
        return intervalTypes[i];
    }

    // Get the tree interval times, types and number of lineages per interval
    protected void updateIntervals() {
        Tree tree = treeInput.get();

        // 1) Get node order by time
        double      [] times = new double[intervalTimes.length];
        SegmentType [] types = new SegmentType[intervalTimes.length];

        Node[] nodes = tree.getNodesAsArray();
        for (int i = 0; i < nodes.length; i++) {
            times[i] = nodes[i].getHeight();
            types[i] = nodes[i].isLeaf() ? SegmentType.SAMPLE : SegmentType.COALESCENT;
        }

        /**** SORTING ****/
        int [] indices = new int[intervalTimes.length];
        HeapSort.sort(times, indices);


        // 2) Traverse nodes from present to TMRCA while reordering and also calculate lineageCounts and intervalWidths
        int lineageCount  = 1;
        intervalTimes[0]  = times[indices[0]];
        intervalTypes[0]  = types[indices[0]];
        intervalWidths[0] = 0;
        lineageCounts[0]  = 0;
        for (int i = 1; i < intervalTimes.length; i++) {

            intervalTimes[i] = times[indices[i]];
            intervalTypes[i] = types[indices[i]];
            intervalWidths[i] = intervalTimes[i]-intervalTimes[i-1];
            lineageCounts[i]  = lineageCount;

            switch (intervalTypes[i]) {
                case SAMPLE:
                    lineageCount++;
                    break;

                case COALESCENT:
                    lineageCount--;
                    break;

                default:
                    break;
            }
        }        
        intervalsKnown = true;

    }
    
    
    /****************************/
    /* Calculation Node methods */
    /****************************/

    @Override
    protected boolean requiresRecalculation() {
        intervalsKnown = false;
        return true;
    }

    @Override
    public void store() {
        //intervalsKnown = false;

        System.arraycopy(intervalTimes, 0, storedIntervalTimes, 0, intervalTimes.length);
        System.arraycopy(intervalTypes, 0, storedIntervalTypes, 0, intervalTypes.length);
        System.arraycopy(intervalWidths, 0, storedIntervalWidths, 0, intervalWidths.length);
        System.arraycopy(lineageCounts, 0, storedLineageCounts, 0, lineageCounts.length);

        super.store();
    }

    @Override
    public void restore() {
        //intervalsKnown = false;

        double [] tmp2 = storedIntervalTimes;
        storedIntervalTimes = intervalTimes;
        intervalTimes = tmp2;

        double [] tmp3 = storedIntervalWidths;
        storedIntervalWidths = intervalWidths;
        intervalWidths = tmp3;

        int [] tmp4 = storedLineageCounts;
        storedLineageCounts = lineageCounts;
        lineageCounts = tmp4;

        SegmentType [] tmp5 = storedIntervalTypes;
        storedIntervalTypes = intervalTypes;
        intervalTypes = tmp5;

        super.restore();
    }

}
