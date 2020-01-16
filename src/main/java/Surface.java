public class Surface {
	IntegralPoint [] pointsOnSurface;
	double detJ;

	public Surface(double ksi1, double eta1, double ksi2, double eta2, double detJ) {
		this.pointsOnSurface = new IntegralPoint[2];
		this.pointsOnSurface[0] = new IntegralPoint(ksi1,eta1,1,1);
		this.pointsOnSurface[1] = new IntegralPoint(ksi2,eta2,1,1);
		this.detJ = detJ;

	}

	public IntegralPoint[] getPointsOnSurface() {
		return pointsOnSurface;
	}
}
