/***
  This file is part of mc3kit.
  
  Copyright (C) 2013 Edward B. Baskerville

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
***/

package mc3kit.util;

import static java.lang.Math.*;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.random.*;
import cern.jet.random.engine.RandomEngine;
import cern.jet.stat.Gamma;

public final class Math {
  public static double LOG_PI = log(java.lang.Math.PI);
  public static double LOG_TWO_PI = log(2 * java.lang.Math.PI);

  public static double logGamma(double x) {
    if(x == 1.0)
      return 0.0;
    return Gamma.logGamma(x);
  }

  public static double logBeta(double x, double y) {
    return logGamma(x) + logGamma(y) - logGamma(x + y);
  }

  public static int[] getRandomPermutation(int size, Uniform uniform) {
    int[] vals = new int[size];
    for(int i = 0; i < size; i++)
      vals[i] = i;
    shuffleInPlace(vals, uniform);
    return vals;
  }

  public static void shuffleInPlace(int[] list, Uniform uniform) {
    for(int i = 0; i < list.length - 1; i++) {
      int tmp = list[i];
      int j = uniform.nextIntFromTo(i, list.length - 1);
      list[i] = list[j];
      list[j] = tmp;
    }
  }

  public static boolean shouldAcceptMetropolisHastings(RandomEngine rng,
      double priorHeatExponent, double likelihoodHeatExponent,
      double oldLogPrior, double oldLogLikelihood, double newLogPrior,
      double newLogLikelihood, double logProposalRatio) {
    boolean accepted = false;
    if(!Double.isInfinite(newLogPrior) && !Double.isNaN(newLogPrior)
        && !Double.isInfinite(newLogLikelihood)
        && !Double.isNaN(newLogLikelihood)) {
      double logPriorRatio = newLogPrior - oldLogPrior;
      double logLikelihoodRatio = newLogLikelihood - oldLogLikelihood;

      if(priorHeatExponent == Double.POSITIVE_INFINITY
          && likelihoodHeatExponent == Double.POSITIVE_INFINITY) {
        if(logPriorRatio + logLikelihoodRatio > 0)
          accepted = true;
      }
      else if(priorHeatExponent == Double.POSITIVE_INFINITY) {
        if(logPriorRatio > 0)
          accepted = true;
      }
      else if(likelihoodHeatExponent == Double.POSITIVE_INFINITY) {
        if(logLikelihoodRatio > 0)
          accepted = true;
      }
      else {
        double logAcceptanceProbability = logProposalRatio + priorHeatExponent
            * logPriorRatio + likelihoodHeatExponent * logLikelihoodRatio;

        if(logAcceptanceProbability >= 0.0
            || log(rng.nextDouble()) <= logAcceptanceProbability) {
          accepted = true;
        }
      }
    }

    return accepted;
  }

  public static double logSumExp(double[] values) {
    double shift = Double.MIN_VALUE;
    for(int i = 0; i < values.length; i++) {
      if(values[i] > shift)
        shift = values[i];
    }

    double sumExpValues = 0;
    for(int i = 0; i < values.length; i++) {
      sumExpValues += exp(values[i] - shift);
    }
    return shift + log(sumExpValues);
  }

  public static double logSumExp(double[] values, double[] coeffs) {
    double shift = Double.MIN_VALUE;
    for(int i = 0; i < values.length; i++) {
      if(values[i] > shift)
        shift = values[i];
    }

    double sumExpValues = 0;
    for(int i = 0; i < values.length; i++) {
      sumExpValues += exp(values[i] - shift) * coeffs[i];
    }
    return shift + log(sumExpValues);
  }

  public static double logSumExp(double[] values, boolean[] negative) {
    double shift = Double.NEGATIVE_INFINITY;
    for(int i = 0; i < values.length; i++) {
      assert !Double.isNaN(values[i]);
      if(Double.isInfinite(values[i])) {
        assert values[i] < 0.0;
        continue;
      }

      if(values[i] > shift) {
        shift = values[i];
      }
    }

    double sumExpValues = 0;
    for(int i = 0; i < values.length; i++) {
      if(Double.isInfinite(values[i]))
        continue;
      sumExpValues += exp(values[i] - shift) * (negative[i] ? -1.0 : 1.0);
    }
    return shift + log(sumExpValues);
  }

  public static double adjustTuningParameter(double param, double measuredRate,
      double targetRate) {
    if(measuredRate == 0) {
      return param / 2;
    }
    else if(measuredRate == 1) {
      return param * 2;
    }

    // Infer target rate using linear interpolation between
    // (lambda = 0, rate = 1) and (lambda = lambda, rate = measuredRate):
    // lambda is updated to be the value on that line where rate = targetRate.
    return param * (1.0 - targetRate) / (1.0 - measuredRate);
  }

  public static double mean(double[] x) {
    double sum = sum(x);
    return sum / x.length;
  }

  public static double center(double[] x) {
    double mean = mean(x);
    for(int i = 0; i < x.length; i++)
      x[i] -= mean;
    return mean;
  }

  public static double sumOfSquares(double[] x) {
    double sumSq = 0.0;
    for(int i = 0; i < x.length; i++)
      sumSq += x[i] * x[i];
    return sumSq;
  }

  public static double norm2(double[] x) {
    return sqrt(sumOfSquares(x));
  }

  public static double norm2(DoubleMatrix1D x) {
    return sqrt(Algebra.DEFAULT.norm2(x));
  }

  public static double normalize(double[] x) {
    double sd = stdDev(x);
    for(int i = 0; i < x.length; i++)
      x[i] /= sd;
    return sd;
  }

  public static double stdDev(double[] x) {
    double mean = mean(x);

    double sumSqErr = 0.0;
    for(int i = 0; i < x.length; i++) {
      double diff = x[i] - mean;
      sumSqErr += diff * diff;
    }

    return sqrt(sumSqErr / x.length);
  }

  public static double stdDev(DoubleMatrix1D x) {
    double mean = mean(x);

    double sumSqErr = 0.0;
    for(int i = 0; i < x.size(); i++) {
      double diff = x.getQuick(i) - mean;
      sumSqErr += diff * diff;
    }

    return sqrt(sumSqErr / x.size());
  }

  public static void divideInPlace(DoubleMatrix1D x, double factor) {
    for(int i = 0; i < x.size(); i++) {
      x.setQuick(i, x.getQuick(i) / factor);
    }
  }

  public static double mean(DoubleMatrix1D x) {
    double sum = 0.0;
    for(int i = 0; i < x.size(); i++)
      sum += x.getQuick(i);
    return sum / x.size();
  }

  public static void standardize(DoubleMatrix1D x) {
    double mean = mean(x);
    double sd = stdDev(x);
    for(int i = 0; i < x.size(); i++) {
      x.setQuick(i, (x.getQuick(i) - mean) / sd);
    }
  }

  public static DoubleMatrix1D subtract(DoubleMatrix1D x1, DoubleMatrix1D x2) {
    assert (x1.size() == x2.size());
    DoubleMatrix1D diff = new DenseDoubleMatrix1D(x1.size());
    for(int i = 0; i < x1.size(); i++)
      diff.setQuick(i, x1.getQuick(i) - x2.getQuick(i));
    return diff;
  }

  // Project x onto p, which is assumed to already be normalized
  public static DoubleMatrix1D project(DoubleMatrix1D x, DoubleMatrix1D p) {
    assert (x.size() == p.size());

    double magnitude = x.zDotProduct(p);
    DoubleMatrix1D xp = new DenseDoubleMatrix1D(x.size());
    for(int i = 0; i < x.size(); i++)
      xp.setQuick(i, p.getQuick(i) * magnitude);

    return xp;
  }

  public static double min(double... x) {
    double min = Double.MIN_VALUE;
    for(double xi : x)
      if(xi < min)
        min = xi;
    return min;
  }

  public static double logisticSigmoid(double x) {
    return 1.0 / (1.0 + exp(-x));
  }

  public static int sum(int... values) {
    int sum = 0;
    for(int value : values)
      sum += value;
    return sum;
  }

  public static double sum(double... values) {
    double sum = 0;
    for(double value : values)
      sum += value;
    return sum;
  }
  
  public static double integrateTrapezoid(double[] x, double[] y) {
    double sum = 0.0;
    for(int i = 0; i < x.length - 1; i++) {
      if(x[i+1] < x[i]) {
        throw new IllegalArgumentException("Integrator assumes x-values are in order smallest to largest.");
      }
      
      sum += (x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2.0;
    }
    return sum;
  }
}
