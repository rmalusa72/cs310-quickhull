class Vector{
	
    int dim;
	float[] values;

	public Vector(int d){
		dim = d;
        values = new float[d];
	}

    public Vector(float[] values){
        dim = values.length;
        this.values = values;
    }

}