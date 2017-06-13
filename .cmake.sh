#!/usr/bin/bash

# Converts the source code of each .cl file insde Library/cl to an string and
# inside an includeable, identical named file in Library

PROJECT_DIR=${PWD}

cd "${PROJECT_DIR}/Library" || exit 1

for CL_SOURCE in ${PROJECT_DIR}/Library/cl/*.cl
do
  CL_BASENAME=${CL_SOURCE##*/}
  CL_FILE="${PROJECT_DIR}/Library/$CL_BASENAME"

  echo "${CL_SOURCE}"

  if [ -f "$CL_FILE" ]
  then
    rm "$CL_FILE"  # Delete existing source file
  fi

  touch "$CL_FILE" # Create new source file

  i=1
  while IFS= read -r LINE
  do
    if (( i < 12 )) # Header part (license, Doxygen file meta)
    then
      printf "%s\n" "$LINE" >> "$CL_FILE"
    elif (( i == 12 )) # #include guard and string name
    then
      CL_MODULE=${CL_BASENAME%.*}
      printf "\n%s\n%s\n\n%s\n%s\n\n%s\n%s\n%s\n%s\n" \
        "#ifndef $(echo "${CL_MODULE}" | tr "[:lower:]" "[:upper:]")_CL" \
        "#define $(echo "${CL_MODULE}" | tr "[:lower:]" "[:upper:]")_CL" \
        "/** @addtogroup ${CL_MODULE/cl/}" \
        " *  @{ */" \
        "/*  @brief The source code to performe operations on the module" \
        " *         @ref ${CL_MODULE/cl/} via OpenCL as a string. */" \
        "static const char ${CL_MODULE}_src[] =" \
        "{" >> "$CL_FILE"
    # Removing include macros and function documentation
  elif [[ $LINE != @(\#include*|/\**| \**) ]]
    then
      # Source code as string #
      LINE=${LINE/\\/\\\\} # Escape backslashes
      LINE=${LINE//\"/\\\"} # Escape quotation marks
      #echo "  \"$LINE\\n\"" >> "$CL_FILE"
      printf "  \"%s\\\n\"\n"  "$LINE" >> "$CL_FILE"
    fi

    let ++i
  done < "$CL_SOURCE"

  # Footer
  printf "%s\n\n%s\n\n%s" \
         "};" \
         "/** @} */" \
         "#endif // $(echo "${CL_BASENAME%.*}" | \
           tr "[:lower:]" "[:upper:]")_CL" >> "$CL_FILE"
done
