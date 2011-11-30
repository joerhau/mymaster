package org.phymod.io;

import org.phymod.Alignment;
import org.phymod.Partition;
import org.phymod.Taxa;

public class PhylipLoader extends AlignmentLoader {
	private FileLoader data;
	private PhylipPartitionLoader part;

	/**
	 * standard contructor
	 * @param l
	 */
	public PhylipLoader(String dataFile, String partFile) {
		super();
		part = new PhylipPartitionLoader(partFile);
		part.parse();
		data = new FileLoader(dataFile);
		data.open();
		
		init();
		alignment = new Alignment(spec);
	}
	
	@Override
	public void init() {
		Partition[] partitions;
		nrPartitions = part.nrPartitions;	
		
		// for each taxa do
		for(int i = 0; i < data.l.size() - 1; i++) {
			partitions = new Partition[part.nrPartitions];
			
			String d = data.line(i + 1).split(" ")[1].replaceAll("\\s+", " ");
			// for each partition do
			for (int j = 0; j < nrPartitions; j++)
				partitions[j] = new Partition(d.substring(part.start[j] - 1, part.end[j]), part.names[j] + "." + part.start[j] + "-" + part.end[j], part.models[j]);
			
			// extract species name
			String name = data.line(i + 1).split(" ")[0].replaceAll("\\s+", " ");
			spec.add(new Taxa(name, partitions));
		}
	}
	
	@Override
	public Alignment getAlignment() {
		return alignment;
	}
}
