import os, glob
from datetime import datetime


class SampleContainer:
    """
    Holds information about the sequencing reads for a particular preprocessing of a sample,
    as well as metadata about the sample's dataset.
    """
    def __init__(self, sample_id, dataset_name, fs_prefix, preprocessing, r1_path, r2_path, last_modified, readable_time):
        self.sample_id = sample_id
        self.dataset_name = dataset_name
        self.fs_prefix = fs_prefix
        self.preprocessing = preprocessing
        self.r1_path = r1_path
        self.r2_path = r2_path
        self.last_modified = last_modified
        self.readable_time = readable_time

    def __str__(self):
        return (f"SampleContainer(sample_id={self.sample_id}, dataset_name={self.dataset_name}, fs_prefix={self.fs_prefix}, "
                f"preprocessing={self.preprocessing}, r1={self.r1_path}, r2={self.r2_path}, "
                f"readable_time={self.readable_time})")

    def __repr__(self):
        return (f"SampleContainer(sample_id={self.sample_id}, dataset_name={self.dataset_name}, fs_prefix={self.fs_prefix}, "
                f"preprocessing={self.preprocessing})")



class Sample:
    """
    Represents a biological sample with its sequencing reads and preprocessing states.
    """
    def __init__(self, sample_id, dataset_path):
        self.sample_id = sample_id
        # Assuming dataset_path is in the form of /path/to/fs_prefix/dataset_name
        self.dataset_path = dataset_path
        self.dataset_name = os.path.basename(dataset_path)  # Extract dataset_name
        self.fs_prefix = os.path.dirname(dataset_path)  # Extract fs_prefix up to the dataset directory
        
        self.containers = self._find_preprocessed_data()

    def _find_preprocessed_data(self):
        containers = []
        dataset_reads_dir = os.path.join(self.dataset_path, 'reads')
        for preproc_dir_name in os.listdir(dataset_reads_dir):
            r1_files = glob.glob(os.path.join(dataset_reads_dir, preproc_dir_name, f"{self.sample_id}_R1*.fastq.gz"))
            r2_files = glob.glob(os.path.join(dataset_reads_dir, preproc_dir_name, f"{self.sample_id}_R2*.fastq.gz"))
            if r1_files and r2_files:
                last_modified_time = os.path.getmtime(r1_files[0])
                readable_time = datetime.fromtimestamp(last_modified_time).strftime("%Y-%b-%d %H:%M:%S")
                container = SampleContainer(self.sample_id, self.dataset_name, self.fs_prefix, preproc_dir_name, r1_files[0], r2_files[0], last_modified_time, readable_time)
                containers.append(container)
        return containers

    def get_latest_container(self):
        return max(self.containers, key=lambda x: x.last_modified, default=None)

    def get_container_by_preprocessing(self, preprocessing):
        return next((x for x in self.containers if x.preprocessing == preprocessing), None)
    
    def __str__(self):
        return f"Sample ID: {self.sample_id}, Dataset Path: {self.dataset_path}, Containers: {len(self.containers)}"
    def __repr__(self):
        return f"Sample ID: {self.sample_id}, Containers: {self.containers}"
