

/**
 * Parameters used for icp initialization
 */
public class ParameterICP
{
        public boolean refine_clamp;         // gleiches wie clip. aber über statistisches Mass, z. B. > 3x mittlerer Abstand, definiert
        public double refine_clamp_par;
        public boolean refine_sd;         // gleiches wie clip. aber über statistisches Mass, z. B. > 3x sd + mittlerer Abstand, definiert
        public double refine_sd_par;
        public boolean refine_clip;          // clipping ab bestimmter Position im Abstandshistogram, 
                                                // z. B. nur 75 % der niedrigeren Abstände berücksichtigt, 
                                                // Rest der Punktepaare verworfen
        public double refine_clip_par;
        
        public boolean refine_sparse;        // nur jeder  xte Punkte wird verwendet
        public double refine_sparse_par;
        public boolean refine_unique;        // klingt als ob jeder Punkt berücksichtigt wird
        public int minimum_valid_points;     // Mindestanzahl gültiger Pixel;
    
	
	public ParameterICP(boolean refine_clamp, double refine_clamp_par, boolean refine_sd, double refine_sd_par, boolean refine_clip, double refine_clip_par, boolean refine_sparse, double refine_sparse_par, boolean refine_unique, int minimum_valid_points)
        {
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
	}
        
    public boolean getRefineClamp(){
        return refine_clamp;
    }
    
    public boolean getRefineSd(){
        return refine_sd;
    }
    
    public boolean getRefineClip(){
        return refine_clip;
    }
    
    public boolean getRefineSparse(){
        return refine_sparse;
    }

    public boolean getRefineUnique(){
        return refine_unique;
    }

    public double getRefineClampPar(){
        return refine_clamp_par;
    }
    
    public double getRefineSdPar(){
        return refine_sd_par;
    }
    
    public double getRefineClipPar(){
        return refine_clip_par;
    }
    
    public int getMinValidPoints(){
        return minimum_valid_points;
    }
    
    public double getRefineSparsePar(){
        return refine_sparse_par;
    }
    
    
   
     public void setRefineClamp(boolean refine_clamp){
        this.refine_clamp = refine_clamp;
    }
     
    public void setRefineSd(boolean refine_sd){
        this.refine_sd = refine_sd;
    }
    
    public void setRefineClip(boolean refine_clip){
        this.refine_clip = refine_clip;
    }
    
    public void setRefineSparse(boolean refine_sparse){
        this.refine_sparse = refine_sparse;
    }

    public void setRefineUnique(boolean refine_unique){
        this.refine_unique = refine_unique;
    }

    public void setRefineClampPar(double refine_clamp_par){
        this.refine_clamp_par = refine_clamp_par;
    }
    
    public void setRefineSdPar(double refine_sd_par){
        this.refine_sd_par = refine_sd_par;
    }
    
    public void setRefineClipPar(double refine_clip_par){
        this.refine_clip_par = refine_clip_par;
    }
    
    public void setRefineSparsePar(double refine_sparse_par){
        this.refine_sparse_par = refine_sparse_par;
    }
    
    public void setMinValidPoints(int minimum_valid_points){
        this.minimum_valid_points = minimum_valid_points;
    }

    @Override
    public String toString()
	{
		return "refine_clamp:\t" + refine_clamp + "\tParameter = " + refine_clamp_par + "\n" +
                        "refine_clip:\t" + refine_clip + "\tParameter  = " + refine_clip_par + "\n" +
                        "refine_sparse:\t" + refine_sparse + "\tParameter = " + refine_sparse_par + "\n" +
                        "refine_unique:\t" + refine_unique;
	}
    
}