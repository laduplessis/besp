package bsp.distributions;

/**
 * Specifies the interval types.
 *
 *
 * Possible other segments:
 *  - Sampled ancestor (like nothing)
 *  - Epoch end time
 *  - Last sample in Epoch
 *  - Polytomy
 *
 */
public enum SegmentType {

    /**
     * Denotes a segment that ends with a coalescent/branching event
     * (In the next segment lineageCount is one smaller)
     */
    COALESCENT("coalescent"),

    /**
     * Denotes a segment that ends with the addition of a new sample
     * (In the next segment lineageCount is one bigger)
     */
    SAMPLE("sample"),

    /**
     * Denotes a segment at the end of which nothing is
     * observed (i.e. the number of lineages is the same in the next interval).
     */
    NOTHING("nothing");

    /**
     * private constructor.
     *
     * @param name the name of the interval type
     */
    private SegmentType(String name) {
        this.name = name;
    }

    @Override
    public String toString() {
            return name;
        }

    private final String name;

}
