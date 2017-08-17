from os.path import join as pjoin
import random

class PATHOLOGYconfig(object): #at some point use pjoin throughout
    def __init__(self):

        #verbosity paramater: if true, will print data about all potentially
        #malignant patches in some zoom level
        self.verbosity = True

        #directory parameters
        self.input_xml = '/Users/doctorfreud/Desktop/serre/Prostate_Bx_Test/Prostate_bx_Test.xml'
        self.svs = '/Users/doctorfreud/Desktop/serre/Prostate_Bx_Test/Prostate_bx_Test.svs'
        self.output_malignant = '/Users/doctorfreud/Desktop/serre/pathology/patches/malignant/'
        self.output_benign = '/Users/doctorfreud/Desktop/serre/pathology/patches/benign/'
        self.home_folder = '/Users/doctorfreud/Desktop/serre/pathology/'

        #patch parameters
        self.patch_size = 254
        self.patch_format = 'JPEG'

        #desired number of benign patches per region per zoom level
        self.number_random_patches = 100
        #desired minimum ratio of tissue/whitespace in random patches
        self.random_whitespace = 0.65 #more tissue than whitespace

        #desired number of benign patches per region per zoom level
        self.number_benign_patches = 100
        #desired minimum ratio of tissue/whitespace in benign patches
        self.benign_whitespace = 0.65
        #maximum allowed malignancy of a 'benign' patch
        self.benign_maximum_malignancy = 1
        #benign w/respect to all annotation together or in separate folders per annotation?
        self.benign_combine_annotations = True

        #annotation parameters
        num_annotation_colors = 100
        annotation_colors = []
        annotation_colors.append('#008000') #1st color: green
        annotation_colors.append('#FF0000') #2nd color: red
        for i in range (2,num_annotation_colors):
            annotation_colors.append("#%06x" % random.randint(0, 0xFFFFFF))
        self.annotation_colors = annotation_colors
