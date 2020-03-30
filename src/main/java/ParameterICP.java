
/**
 * Parameters used for icp initialization
 */
class ParameterICP {
	private boolean refine_clamp; // gleiches wie clip. aber über statistisches Mass, z. B. > 3x mittlerer
									// Abstand, definiert
	private double refine_clamp_par;
	private boolean refine_sd; // gleiches wie clip. aber über statistisches Mass, z. B. > 3x sd + mittlerer
								// Abstand, definiert
	private double refine_sd_par;
	private boolean refine_clip; // clipping ab bestimmter Position im Abstandshistogram,
									// z. B. nur 75 % der niedrigeren Abstände berücksichtigt,
									// Rest der Punktepaare verworfen
	private double refine_clip_par;

	private boolean refine_sparse; // nur jeder xte Punkte wird verwendet
	private double refine_sparse_par;
	private boolean refine_unique; // klingt als ob jeder Punkt berücksichtigt wird
	private int minimum_valid_points; // Mindestanzahl gültiger Pixel;

	private boolean refine_a_priori_landmark;
	private boolean refine_a_priori_diff_slider;
	private boolean use_landmark_mask;

	public ParameterICP(boolean refine_clamp, double refine_clamp_par, boolean refine_sd, double refine_sd_par,
						boolean refine_clip, double refine_clip_par, boolean refine_sparse, double refine_sparse_par,
						boolean refine_unique, int minimum_valid_points, boolean refine_a_priori_landmark,
						boolean refine_a_priori_diff_slider, boolean use_landmark_mask) {
		this.refine_clamp = refine_clamp;
		this.refine_clamp_par = refine_clamp_par;
		this.refine_sd = refine_sd;
		this.refine_sd_par = refine_sd_par;
		this.refine_clip = refine_clip;
		this.refine_clip_par = refine_clip_par;
		this.refine_sparse = refine_sparse;
		this.refine_sparse_par = refine_sparse_par;
		this.refine_unique = refine_unique;
		this.minimum_valid_points = minimum_valid_points;
		this.refine_a_priori_landmark = refine_a_priori_landmark;
		this.refine_a_priori_diff_slider = refine_a_priori_diff_slider;
		this.use_landmark_mask = use_landmark_mask;
	}

	public boolean getRefineClamp() {
		return refine_clamp;
	}

	public void setRefineClamp(boolean refine_clamp) {
		this.refine_clamp = refine_clamp;
	}

	public boolean getRefineSd() {
		return refine_sd;
	}

	public void setRefineSd(boolean refine_sd) {
		this.refine_sd = refine_sd;
	}

	public boolean getRefineClip() {
		return refine_clip;
	}

	public void setRefineClip(boolean refine_clip) {
		this.refine_clip = refine_clip;
	}

	public boolean getRefineSparse() {
		return refine_sparse;
	}

	public void setRefineSparse(boolean refine_sparse) {
		this.refine_sparse = refine_sparse;
	}

	public boolean getRefineUnique() {
		return refine_unique;
	}

	public void setRefineUnique(boolean refine_unique) {
		this.refine_unique = refine_unique;
	}

	public double getRefineClampPar() {
		return refine_clamp_par;
	}

	public void setRefineClampPar(double refine_clamp_par) {
		this.refine_clamp_par = refine_clamp_par;
	}

	public double getRefineSdPar() {
		return refine_sd_par;
	}

	public void setRefineSdPar(double refine_sd_par) {
		this.refine_sd_par = refine_sd_par;
	}

	public double getRefineClipPar() {
		return refine_clip_par;
	}

	public void setRefineClipPar(double refine_clip_par) {
		this.refine_clip_par = refine_clip_par;
	}

	public int getMinValidPoints() {
		return minimum_valid_points;
	}

	public void setMinValidPoints(int minimum_valid_points) {
		this.minimum_valid_points = minimum_valid_points;
	}

	public double getRefineSparsePar() {
		return refine_sparse_par;
	}

	public void setRefineSparsePar(double refine_sparse_par) {
		this.refine_sparse_par = refine_sparse_par;
	}

	public boolean getRefineAPrioriLandmark() {
		return refine_a_priori_landmark;
	}

	public void setRefineAPrioriLandmark(boolean refine_a_priori_landmark) {
		this.refine_a_priori_landmark = refine_a_priori_landmark;
	}

	public boolean getRefineAPrioriDiffSlider() {
		return refine_a_priori_diff_slider;
	}

	public void setRefineAPrioriDiffSlider(boolean refine_a_priori_diff_slider) {
		this.refine_a_priori_diff_slider = refine_a_priori_diff_slider;
	}

	public boolean getUseLandmarkMask() {
		return use_landmark_mask;
	}

	public void setUseLandmarkMask(boolean use_landmark_mask) {
		this.use_landmark_mask = use_landmark_mask;
	}

	@Override
	public String toString() {
		return "refine_clamp:\t" + refine_clamp + "\tParameter = " + refine_clamp_par + "\n" + "refine_clip:\t"
				+ refine_clip + "\tParameter  = " + refine_clip_par + "\n" + "refine_sparse:\t" + refine_sparse
				+ "\tParameter = " + refine_sparse_par + "\n" + "refine_unique:\t" + refine_unique
				+ "refine_a_priori_landmark:\t" + refine_a_priori_landmark + "\n" + "refine_a_priori_diff_slider:\t"
				+ refine_a_priori_diff_slider + "\n"+ "use_landmark_mask:\t" + use_landmark_mask;
	}

}