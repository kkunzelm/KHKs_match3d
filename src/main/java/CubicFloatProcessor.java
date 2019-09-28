import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
/** Version of a FloatProcessor with constrained cubic interpolation */
// Bob Dougherty 6/19/2008. Modified from FloatProcesor by Wayne Rasband.
// http://www.optinav.com/CubicFloatResizeRotate.htm
class CubicFloatProcessor extends FloatProcessor {
	private float[] pixels;
	/**
	 * Creates a blank CubicFloatProcessor using the default grayscale LUT that
	 * displays zero as black. Call invertLut() to display zero as white.
	 */
	private CubicFloatProcessor(int width, int height) {
		super(width, height);
		pixels = (float[]) getPixels();
	}

	/**
	 * Creates a CubicFloatProcessor from an int array using the default grayscale
	 * LUT.
	 */
	public CubicFloatProcessor(int width, int height, int[] pixelsInt) {
		super(width, height, pixelsInt);
		pixels = (float[]) getPixels();
	}

	/**
	 * Creates a CubicFloatProcessor from a double array using the default grayscale
	 * LUT.
	 */
	public CubicFloatProcessor(int width, int height, double[] pixelsDouble) {
		super(width, height, pixelsDouble);
		pixels = (float[]) getPixels();
	}

	/** Creates a FloatProcessor from a float[][] array using the default LUT. */
	public CubicFloatProcessor(float[][] array) {
		super(array);
		pixels = (float[]) getPixels();
	}

	/** Creates a CubicFloatProcessor from an int[][] array. */
	public CubicFloatProcessor(int[][] array) {
		super(array);
		pixels = (float[]) getPixels();
	}
	public void setPixels(Object pixelsNew) {
		super.setPixels(pixelsNew);
		pixels = (float[]) getPixels();
	}
	/**
	 * Creates a new CubicFloatProcessor containing a scaled copy of this image or
	 * selection.
	 */

	public ImageProcessor resize(int dstWidth, int dstHeight) {
		float[] pixels = (float[]) getPixels();
		double srcCenterX = roiX + roiWidth / 2.0;
		double srcCenterY = roiY + roiHeight / 2.0;
		double dstCenterX = dstWidth / 2.0;
		double dstCenterY = dstHeight / 2.0;
		double xScale = (double) dstWidth / roiWidth;
		double yScale = (double) dstHeight / roiHeight;
		if (interpolate) {
			dstCenterX += xScale / 2.0;
			dstCenterY += yScale / 2.0;
		}
		CubicFloatProcessor ip2 = new CubicFloatProcessor(dstWidth, dstHeight);
		float[] pixels2 = (float[]) ip2.getPixels();
		double xs, ys;
		double xlimit = width - 1.0, xlimit2 = width - 1.001;
		double ylimit = height - 1.0, ylimit2 = height - 1.001;
		int index1, index2;
		for (int y = 0; y <= dstHeight - 1; y++) {
			ys = (y - dstCenterY) / yScale + srcCenterY;
			if (interpolate) {
				if (ys < 0.0)
					ys = 0.0;
				if (ys >= ylimit)
					ys = ylimit2;
			}
			index1 = width * (int) ys;
			index2 = y * dstWidth;
			for (int x = 0; x <= dstWidth - 1; x++) {
				xs = (x - dstCenterX) / xScale + srcCenterX;
				if (interpolate) {
					if (xs < 0.0)
						xs = 0.0;
					if (xs >= xlimit)
						xs = xlimit2;
					pixels2[index2++] = (float) getInterpolatedPixel(xs, ys, pixels);
				} else
					pixels2[index2++] = pixels[index1 + (int) xs];
			}
		}
		return ip2;
	}
	/**
	 * Uses cubic interpolation to find the pixel value at real coordinates (x,y).
	 */
	public double getInterpolatedPixel(double x, double y) {
		if (x < 0.0)
			x = 0.0;
		if (x >= width - 1.0)
			x = width - 1.001;
		if (y < 0.0)
			y = 0.0;
		if (y >= height - 1.0)
			y = height - 1.001;
		return getInterpolatedPixel(x, y, pixels);
	}
	/**
	 * Rotates the image or ROI 'angle' degrees clockwise.
	 * 
	 * @see ImageProcessor#setInterpolate
	 */
	public void rotate(double angle) {
		float[] pixels2 = (float[]) getPixelsCopy();
		double centerX = roiX + (roiWidth - 1) / 2.0;
		double centerY = roiY + (roiHeight - 1) / 2.0;
		int xMax = roiX + this.roiWidth - 1;

		double angleRadians = -angle / (180.0 / Math.PI);
		double ca = Math.cos(angleRadians);
		double sa = Math.sin(angleRadians);
		double tmp1 = centerY * sa - centerX * ca;
		double tmp2 = -centerX * sa - centerY * ca;
		double tmp3, tmp4, xs, ys;
		int index, ixs, iys;
		double dwidth = width, dheight = height;
		double xlimit = width - 1.0, xlimit2 = width - 1.001;
		double ylimit = height - 1.0, ylimit2 = height - 1.001;

		for (int y = roiY; y < (roiY + roiHeight); y++) {
			index = y * width + roiX;
			tmp3 = tmp1 - y * sa + centerX;
			tmp4 = tmp2 + y * ca + centerY;
			for (int x = roiX; x <= xMax; x++) {
				xs = x * ca + tmp3;
				ys = x * sa + tmp4;
				if ((xs >= -0.01) && (xs < dwidth) && (ys >= -0.01) && (ys < dheight)) {
					if (interpolate) {
						if (xs < 0.0)
							xs = 0.0;
						if (xs >= xlimit)
							xs = xlimit2;
						if (ys < 0.0)
							ys = 0.0;
						if (ys >= ylimit)
							ys = ylimit2;
						pixels[index++] = (float) getInterpolatedPixel(xs, ys, pixels2);
					} else {
						ixs = (int) (xs + 0.5);
						iys = (int) (ys + 0.5);
						if (ixs >= width)
							ixs = width - 1;
						if (iys >= height)
							iys = height - 1;
						pixels[index++] = pixels2[width * iys + ixs];
					}
				} else
					pixels[index++] = 0;
			}
			if (y % 20 == 0)
				showProgress((double) (y - roiY) / roiHeight);
		}
		showProgress(1.0);
	}
	/**
	 * Scales the image or selection using the specified scale factors.
	 * 
	 * @see ImageProcessor#setInterpolate
	 */
	public void scale(double xScale, double yScale) {
		double xCenter = roiX + roiWidth / 2.0;
		double yCenter = roiY + roiHeight / 2.0;
		int xmin, xmax, ymin, ymax;
		double min = getMin();

		if ((xScale > 1.0) && (yScale > 1.0)) {
			// expand roi
			xmin = (int) (xCenter - (xCenter - roiX) * xScale);
			if (xmin < 0)
				xmin = 0;
			xmax = xmin + (int) (roiWidth * xScale) - 1;
			if (xmax >= width)
				xmax = width - 1;
			ymin = (int) (yCenter - (yCenter - roiY) * yScale);
			if (ymin < 0)
				ymin = 0;
			ymax = ymin + (int) (roiHeight * yScale) - 1;
			if (ymax >= height)
				ymax = height - 1;
		} else {
			xmin = roiX;
			xmax = roiX + roiWidth - 1;
			ymin = roiY;
			ymax = roiY + roiHeight - 1;
		}
		float[] pixels2 = (float[]) getPixelsCopy();
		boolean checkCoordinates = (xScale < 1.0) || (yScale < 1.0);
		int index1, index2, xsi, ysi;
		double ys, xs;
		double xlimit = width - 1.0, xlimit2 = width - 1.001;
		double ylimit = height - 1.0, ylimit2 = height - 1.001;
		for (int y = ymin; y <= ymax; y++) {
			ys = (y - yCenter) / yScale + yCenter;
			ysi = (int) ys;
			if (ys < 0.0)
				ys = 0.0;
			if (ys >= ylimit)
				ys = ylimit2;
			index1 = y * width + xmin;
			index2 = width * (int) ys;
			for (int x = xmin; x <= xmax; x++) {
				xs = (x - xCenter) / xScale + xCenter;
				xsi = (int) xs;
				if (checkCoordinates && ((xsi < xmin) || (xsi > xmax) || (ysi < ymin) || (ysi > ymax)))
					pixels[index1++] = (float) min;
				else {
					if (interpolate) {
						if (xs < 0.0)
							xs = 0.0;
						if (xs >= xlimit)
							xs = xlimit2;
						pixels[index1++] = (float) getInterpolatedPixel(xs, ys, pixels2);
					} else
						pixels[index1++] = pixels2[index2 + xsi];
				}
			}
			if (y % 20 == 0)
				showProgress((double) (y - ymin) / height);
		}
		showProgress(1.0);
	}
	/**
	 * Uses bicubic interpolation to find the pixel value at real coordinates (x,y).
	 */
	private double getInterpolatedPixel(double x, double y, float[] pixels) {
		int xBase = (int) x;
		int yBase = (int) y;
		double xFraction = x - xBase;
		double yFraction = y - yBase;
		int offset = yBase * width + xBase;
		double ll = pixels[offset];
		double lr = pixels[offset + 1];
		double ur = pixels[offset + width + 1];
		double ul = pixels[offset + width];
		/*
		 * UpperLeft, LowerRight, etc. The base point is ll xFraction and yFraction are
		 * >= 0 and < 1 It could happen that the two-index points, like l_ll, are off
		 * the image. In this case they will be extrapolated from the inside points.
		 * 
		 * x->
		 * 
		 * ll_ll ll_l ll_r ll_rr y l_ll ll lr l_rr | u_ll ul ur u_rr v uu_ll uu_l uu_r
		 * uu_rr
		 * 
		 */
		// Initially extrapolate everybody to make sure something
		// happens
		double ll_ll = 2 * ll - ur;
		double ll_l = 2 * ll - ul;
		double ll_r = 2 * lr - ur;
		double ll_rr = 2 * lr - ul;
		double l_rr = 2 * lr - ll;
		double u_rr = 2 * ur - ul;
		double uu_rr = 2 * ur - ll;
		double uu_r = 2 * ur - lr;
		double uu_l = 2 * ul - ll;
		double uu_ll = 2 * ul - lr;
		double u_ll = 2 * ul - ur;
		double l_ll = 2 * ll - lr;
		// Look up real pixels if they are in the image.
		// Top row of matrix
		if (yBase > 0) {
			ll_l = pixels[offset - width];
			ll_r = pixels[offset + 1 - width];
			if (xBase > 0) {
				ll_ll = pixels[offset - width - 1];
			}
			if (xBase < (width - 2)) {
				ll_rr = pixels[offset + 2 - width];
			}
		}
		// Bottom row of matrix
		if (yBase < (height - 2)) {
			uu_l = pixels[offset + width + width];
			uu_r = pixels[offset + width + 1 + width];
			if (xBase > 0) {
				uu_ll = pixels[offset + width + width - 1];
			}
			if (xBase < (width - 2)) {
				uu_rr = pixels[offset + width + width + 2];
			}
		}
		// left column
		if (xBase > 0) {
			l_ll = pixels[offset - 1];
			u_ll = pixels[offset + width - 1];
		}
		// right column
		if (xBase < (width - 2)) {
			l_rr = pixels[offset + 2];
			u_rr = pixels[offset + width + 2];
		}
		// Interpolate each of the rows
		double lli = cubic(xFraction, ll_ll, ll_l, ll_r, ll_rr);
		double li = cubic(xFraction, l_ll, ll, lr, l_rr);
		double ui = cubic(xFraction, u_ll, ul, ur, u_rr);
		double uui = cubic(xFraction, uu_ll, uu_l, uu_r, uu_rr);
		// Interpolate in the vertical direction
		return cubic(yFraction, lli, li, ui, uui);
	}
	// Four point cubic interpolation, where it is assumed that the evaluation
	// point is within the center interval of the three intervals. Overshoot
	// is controlled by constraining the result to not differ too much from
	// linear interpolation using just the center points. The amount of
	// allowed difference is determined by performing linear extrapolation
	// from the outside intervals.
	private double cubic(double t, double fm1, double f0, double f1, double f2) {
		double dm1 = fm1 - f0;
		double d1 = f1 - f0;
		double d2 = f2 - f0;
		double b = (d1 + dm1) / 2;
		double d1_minus_b = d1 - b;
		double d2_minus_4b = d2 - 4 * b;
		double c = (d2_minus_4b - 2 * d1_minus_b) / 6;
		double a = d1_minus_b - c;
		double cubic = f0 + t * (a + t * (b + t * c));
		double linear = linInterp(t, f0, f1);
		double extrapLeft = linInterp(t + 1, fm1, f0);
		double extrapRight = linInterp(t - 1, f1, f2);
		// Overshoot control idea: the interpolated answer must be
		// between the linear interpolation result and whichever
		// linear extrapolation result is closer to the linear
		// interpolation result. If the cubic result is not in this
		// interval, then the answer is either the linear or the closer-
		// to-linear extrapolated.
		double closerExtrap = extrapLeft;
		if (Math.abs(extrapRight - linear) < Math.abs(extrapLeft - linear)) {
			closerExtrap = extrapRight;
		}
		double upperBound = closerExtrap;
		double lowerBound = linear;
		if (upperBound < lowerBound) {
			upperBound = linear;
			lowerBound = closerExtrap;
		}
		if (cubic > upperBound) {
			return upperBound;
		} else if (cubic < lowerBound) {
			return lowerBound;
		}
		return cubic;
	}
	private double linInterp(double t, double f0, double f1) {
		return (1 - t) * f0 + t * f1;
	}
}
