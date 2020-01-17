// 
// Decompiled by Procyon v0.5.36
// 

package multifilediversity;

import org.apache.commons.math3.special.Gamma;

public class GammaValuesPoly
{
    double[] alpha;
    double sumAlpha;
    double[] lGammaAlpha;
    double lGammaSumAlpha;
    double[][] lGammaAlphaPlusN;
    double[] lGammaSumAlphaPlusN;
    double[] lGammaN;
    
    GammaValuesPoly(final double[] alpha) {
        this.lGammaAlpha = new double[4];
        this.lGammaAlphaPlusN = new double[4][1000];
        this.lGammaSumAlphaPlusN = new double[1000];
        this.lGammaN = new double[1000];
        this.alpha = alpha;
        this.sumAlpha = alpha[0] + alpha[1] + alpha[2] + alpha[3];
        this.lGammaSumAlpha = Gamma.logGamma(this.sumAlpha);
        for (int iAlpha = 0; iAlpha < 4; ++iAlpha) {
            this.lGammaAlpha[iAlpha] = Gamma.logGamma(alpha[iAlpha]);
        }
        for (int n = 0; n < 1000; ++n) {
            this.lGammaSumAlphaPlusN[n] = Gamma.logGamma(this.sumAlpha + n);
            this.lGammaAlphaPlusN[0][n] = Gamma.logGamma(alpha[0] + n);
            this.lGammaAlphaPlusN[1][n] = Gamma.logGamma(alpha[1] + n);
            this.lGammaAlphaPlusN[2][n] = Gamma.logGamma(alpha[2] + n);
            this.lGammaAlphaPlusN[3][n] = Gamma.logGamma(alpha[3] + n);
            this.lGammaN[n] = Gamma.logGamma((double)n);
        }
    }
    
    double getLGammaAlphaPlusN(final int i, final int n) {
        if (n < 1000) {
            return this.lGammaAlphaPlusN[i][n];
        }
        return Gamma.logGamma(this.alpha[i] + n);
    }
    
    double getLGammaSumAlphaPlusN(final int n) {
        if (n < 1000) {
            return this.lGammaSumAlphaPlusN[n];
        }
        return Gamma.logGamma(this.sumAlpha + n);
    }
    
    double getLGammaN(final int n) {
        if (n < 1000) {
            return this.lGammaN[n];
        }
        return Gamma.logGamma((double)n);
    }
}
