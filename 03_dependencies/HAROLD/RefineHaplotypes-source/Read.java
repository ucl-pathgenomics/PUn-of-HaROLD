package refineHaplotypes;

import java.util.Iterator;
import java.util.ArrayList;
import java.util.Arrays;
import htsjdk.samtools.SAMRecord;
import java.util.HashMap;

public class Read
{
    private int iStart;
    private int iLength;
    private int[] sequence;
    private byte[] baseQual;
    private boolean[] siteExists;
    private int[] limits;
    private boolean negativeStrand;
    private HashMap<Integer, Integer> sigSiteHash;
    private int nSigSites;
    private int nCopies;
    private String sigSiteTag;
    
    public Read(final SAMRecord samRecord) {
        this.limits = new int[2];
        this.sigSiteHash = new HashMap<Integer, Integer>();
        this.nSigSites = 0;
        this.nCopies = 1;
        this.sigSiteTag = "";
        this.iStart = samRecord.getAlignmentStart();
        this.iLength = samRecord.getAlignmentEnd() - samRecord.getAlignmentStart();
        this.negativeStrand = samRecord.getReadNegativeStrandFlag();
        this.limits[0] = this.iStart;
        this.limits[1] = this.iStart + this.iLength;
        this.sequence = new int[this.iLength];
        this.baseQual = new byte[this.iLength];
        this.siteExists = new boolean[this.iLength];
        Arrays.fill(this.sequence, 4);
        Arrays.fill(this.baseQual, (byte)0);
        Arrays.fill(this.siteExists, false);
        for (int iSite = this.iStart; iSite < this.iStart + this.iLength; ++iSite) {
            final int readPosition = samRecord.getReadPositionAtReferencePosition(iSite) - 1;
            if (readPosition >= 0 && (samRecord.getBaseQualities().length == 0 || samRecord.getBaseQualities()[readPosition] > RefineHaplotypes.minBaseQual)) {
                this.siteExists[iSite - this.iStart] = true;
                this.sequence[iSite - this.iStart] = RefineHaplotypes.getDNA(samRecord.getReadString().charAt(readPosition));
            }
        }
    }
    
    public HashMap<Integer, Integer> getSigSiteHash() {
        return this.sigSiteHash;
    }
    
    public void setSignificantSites(final ArrayList<Integer> sigSiteList) {
        for (final int iSite : sigSiteList) {
            final int lSite = iSite - this.iStart;
            if (iSite >= this.iStart && iSite < this.iStart + this.iLength && this.siteExists[lSite] && this.sequence[lSite] < 4) {
                this.sigSiteHash.put(iSite, this.sequence[lSite]);
                if (this.sigSiteTag.length() > 0) {
                    this.sigSiteTag = invokedynamic(makeConcatWithConstants:(Ljava/lang/String;II)Ljava/lang/String;, this.sigSiteTag, iSite, this.sequence[lSite]);
                }
                else {
                    this.sigSiteTag = invokedynamic(makeConcatWithConstants:(Ljava/lang/String;II)Ljava/lang/String;, this.sigSiteTag, iSite, this.sequence[lSite]);
                }
            }
        }
        this.nSigSites = this.sigSiteHash.size();
    }
    
    public String getSigSiteTag() {
        return this.sigSiteTag;
    }
    
    public int getNCopies() {
        return this.nCopies;
    }
    
    public void addCopies() {
        ++this.nCopies;
    }
    
    public int[] getSequence() {
        return this.sequence;
    }
    
    public boolean[] getSiteExists() {
        return this.siteExists;
    }
    
    public int[] getLimits() {
        return this.limits;
    }
    
    public int getBase(final int iSite) {
        return this.sequence[iSite - this.iStart];
    }
    
    public boolean getNegativeStrand() {
        return this.negativeStrand;
    }
} 