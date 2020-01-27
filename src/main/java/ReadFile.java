import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

class ReadFile {
	double initialTemperature;
	double simulationTime;
	double simulationStepTime;
	double ambientTemperature;
	double alpha;
	double height;
	double width;
	int nodesPerHeight;
	int nodesPerWidth;
	double specificHeat;
	double conductivity;
	double density;

	ReadFile(String path) throws IOException {
		List lines = readFileIntoList(path);
		if(lines.size()==12) {
			this.initialTemperature = Double.parseDouble(String.valueOf(lines.get(0)));
			this.simulationTime = Double.parseDouble(String.valueOf(lines.get(1)));
			this.simulationStepTime = Double.parseDouble(String.valueOf(lines.get(2)));
			this.ambientTemperature = Double.parseDouble(String.valueOf(lines.get(3)));
			this.alpha = Double.parseDouble(String.valueOf(lines.get(4)));
			this.height = Double.parseDouble(String.valueOf(lines.get(5)));
			this.width = Double.parseDouble(String.valueOf(lines.get(6)));
			this.nodesPerHeight = Integer.parseInt(String.valueOf(lines.get(7)));
			this.nodesPerWidth = Integer.parseInt(String.valueOf(lines.get(8)));
			this.specificHeat = Double.parseDouble(String.valueOf(lines.get(9)));
			this.conductivity = Double.parseDouble(String.valueOf(lines.get(10)));
			this.density = Double.parseDouble(String.valueOf(lines.get(11)));
		}
		else throw new IllegalArgumentException();
	}

	private static List <String> readFileIntoList(String path) throws IOException {
		File file = new File(path);
		BufferedReader bufferedReader = new BufferedReader(new FileReader(file));
		List <String> lines = new ArrayList<>();
		String line;
		while ((line = bufferedReader.readLine())!=null)
			lines.add(line);
		return lines;
	}

}
