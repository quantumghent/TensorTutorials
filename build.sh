#!/bin/bash
echo "Building the project..."

usage() { echo "Usage: $0 [-h] [-d]" 1>&2; exit 1; }

while getopts 'h:d:' arg; do
    case "${arg}" in
        h)
            usage
            ;;
        d)
            echo "Installing dependencies:"

            echo "Installing Python dependencies"
            pip install -r lectures/requirements.txt

            echo "Installing LaTeX dependencies"
            sudo apt-get -qq update
            sudo apt-get install -y     \
                texlive-latex-recommended \
                texlive-latex-extra       \
                texlive-fonts-recommended \
                texlive-fonts-extra       \
                texlive-xetex             \
                latexmk                   \
                xindy                     \
                texlive-luatex            \
                dvipng                    \
                ghostscript               \
                cm-super

            echo "Install IJulia and Setup Project"
            julia -e 'using Pkg; Pkg.add("IJulia");'
            julia --project=lectures --threads auto -e 'using Pkg; Pkg.instantiate();'
            ;;
        *)
            usage
            ;;
    esac
done

echo "Building the lectures"
jupyter-book build lectures/