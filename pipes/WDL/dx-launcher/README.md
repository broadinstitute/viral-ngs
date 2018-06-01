# dx-launcher for viral-ngs demux

This applet is meant to be invoked automatically by [dx-streaming-upload](https://github.com/dnanexus-rnd/dx-streaming-upload) when it completes upload of a sequencing run. The applet:

1. Reads the "upload sentinel record" to get the run metadata and identify the uploaded tarball(s)
2. If the run has been uploaded in multiple tarballs over time, consolidates them into one.
3. Launches the demux_only/demux_plus WDL workflow for each lane of the run.
