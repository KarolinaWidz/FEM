class Surface {
	private IntegralPoint [] pointsOnSurface;
	double jacobian1Ddet;

	Surface(double ksi1, double eta1, double ksi2, double eta2, double jacobian1Ddet) {
		this.pointsOnSurface = new IntegralPoint[2];
		this.pointsOnSurface[0] = new IntegralPoint(ksi1,eta1,1,1);
		this.pointsOnSurface[1] = new IntegralPoint(ksi2,eta2,1,1);
		this.jacobian1Ddet = jacobian1Ddet;

	}

	public IntegralPoint[] getPointsOnSurface() {
		return pointsOnSurface;
	}
}
