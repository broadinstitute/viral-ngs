# dx-launcher for viral-ngs demux

This DNAnexus applet is meant to be invoked automatically by [dx-streaming-upload](https://github.com/dnanexus-rnd/dx-streaming-upload) when it completes upload of a sequencing run. The applet:

1. Reads the "upload sentinel record" to get the run metadata and identify the uploaded tarball(s)
2. If the run has been uploaded in multiple tarballs over time, consolidates them into one
3. Launches the demux_only/demux_plus WDL workflow for each lane of the run

See [demux_launcher.yml](https://github.com/broadinstitute/viral-ngs/blob/master/pipes/WDL/dx-launcher/demux_launcher.yml) for detailed input documentation. The repository's Travis scripts build and test this applet ([travis/build-dx.sh](https://github.com/broadinstitute/viral-ngs/blob/master/travis/build-dx.sh) [travis/tests-dx.sh](https://github.com/broadinstitute/viral-ngs/blob/mlin-dx-launcher/travis/tests-dx.sh)).

By default, the demux workflows are launched as children of the dx-launcher job. This means that the whole job tree will fail if any one of them does. Optionally, each lane demux can be launched as an independent top-level job, to succeed or fail independently. This requires a DNAnexus API key (with at least Contributor access to the project), provided in a text file as an input to the job. To protect this key, the file should be isolated in a separate DNAnexus project, accessible only to the user account that will invoke dx-launcher (and necessary administrators). The applet can then be launched with e.g. `-i api_token=token_project:/token.txt`

As a side-effect of this workaround, the destination folder for the outputs must be specified using the `folder` input to the applet (`-i folder=/xxx`), instead of the system's folder setting (`--folder /xxx`).
