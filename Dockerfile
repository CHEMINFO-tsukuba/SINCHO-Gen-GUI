# syntax=docker/dockerfile:1

ARG SINCHO_BASE_IMAGE=sincho-gen:0.1.0
FROM ${SINCHO_BASE_IMAGE}

ARG SINCHO_VERSION=0.1.0
LABEL org.opencontainers.image.title="SINCHO-Gen-GUI" \
      org.opencontainers.image.version="${SINCHO_VERSION}"

COPY . /SINCHO-Gen-GUI
RUN python3 -m pip install /SINCHO-Gen-GUI && \
    python3 -m pip install pyside2 stmol && \
    python3 -m pip install --no-deps protobuf==3.20.3 && \
    python3 -c "import google.protobuf, streamlit, pandas, tensorflow; assert google.protobuf.__version__ == '3.20.3'; assert pandas.__version__.startswith('1.')"
    # python3 -c "import streamlit, pandas; assert pandas.__version__.startswith('1.')"

WORKDIR /work
EXPOSE 8501
CMD ["sincho-gen-gui"]
