from Biomarker import Biomarker

datapath = "example_data/Nexus.txt"
# datapath = "example_data/Nihonkoden.txt"
# datapath = "example_data/Biolog.csv"

Bio = Biomarker(datapath, DeviceName="Nexus")
# Bio = Biomarker(datapath, DeviceName="Nihonkoden")
# Bio = Biomarker(datapath, DeviceName="Biolog")

Bio.showGraph()

# print(dir(Bio))
