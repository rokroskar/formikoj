from formikoj import SeismicWaveformModeler

swm = SeismicWaveformModeler('.')

swm.load('mesh')
swm.load('config')
swm.create('waveforms')
