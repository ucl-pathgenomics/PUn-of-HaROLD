// 
// Decompiled by Procyon v0.5.36
// 

package multifilediversity;

import org.apache.commons.math3.special.Gamma;

public class GammaValuesErrors
{
    double alpha0;
    double alphaE;
    double lGammaAlpha0PlusAlphaE;
    double lGammaAlphaE;
    double lGammaAlpha0PlusFourAlphaE;
    double[] lGammaAlpha0PlusAlphaEPlusN;
    double[] lGammaAlphaEPlusN;
    double[] lGammaN;
    double[] lGammaAlpha0PlusFourAlphaEPlusN;
    
    GammaValuesErrors(final double alpha0, final double alphaE) {
        this.lGammaAlpha0PlusAlphaEPlusN = new double[1000];
        this.lGammaAlphaEPlusN = new double[1000];
        this.lGammaN = new double[1000];
        this.lGammaAlpha0PlusFourAlphaEPlusN = new double[1000];
        this.alpha0 = alpha0;
        this.alphaE = alphaE;
        this.lGammaAlphaE = Gamma.logGamma(alphaE);
        this.lGammaAlpha0PlusAlphaE = Gamma.logGamma(alpha0 + alphaE);
        this.lGammaAlpha0PlusFourAlphaE = Gamma.logGamma(alpha0 + 4.0 * alphaE);
        for (int n = 0; n < 1000; ++n) {
            this.lGammaAlphaEPlusN[n] = Gamma.logGamma(alphaE + n);
            this.lGammaAlpha0PlusAlphaEPlusN[n] = Gamma.logGamma(alpha0 + alphaE + n);
            this.lGammaAlpha0PlusFourAlphaEPlusN[n] = Gamma.logGamma(alpha0 + 4.0 * alphaE + n);
            this.lGammaN[n] = Gamma.logGamma((double)n);
        }
    }
    
    double getLGammaAlpha0PlusFourAlphaEPlusN(final int n) {
        if (n < 1000) {
            return this.lGammaAlpha0PlusFourAlphaEPlusN[n];
        }
        return Gamma.logGamma(this.alpha0 + 4.0 * this.alphaE + n);
    }
    
    double getLGammaN(final int n) {
        if (n < 1000) {
            return this.lGammaN[n];
        }
        return Gamma.logGamma((double)n);
    }
    
    double getLGammaAlpha0PlusAlphaEPlusN(final int n) {
        if (n < 1000) {
            return this.lGammaAlpha0PlusAlphaEPlusN[n];
        }
        return Gamma.logGamma(this.alpha0 + this.alphaE + n);
    }
    
    double getLGammaAlphaEPlusN(final int n) {
        if (n < 1000) {
            return this.lGammaAlphaEPlusN[n];
        }
        return Gamma.logGamma(this.alphaE + n);
    }
}
