#!/bin/bash

# Text color codes
RED="\033[0;31m"
GREEN="\033[0;32m"
YELLOW="\033[0;33m"
BLUE="\033[0;34m"
RESET="\033[0m"  # Reset text formatting

# Usage: colored_echo "color_code" "message"
colored_echo() {
  local color="$1"
  local message="$2"
  echo -e "${color}${message}${RESET}"
}

#colored_echo "$RED" "This is a red message."
#colored_echo "$GREEN" "This is a green message."
#colored_echo "$YELLOW" "This is a yellow message."
#colored_echo "$BLUE" "This is a blue message."

get_workdir() {
  local file="$1"
  value=$(jq '.WorkDir' "$file")
  if [[ "$value" == "null" ]]; then
	  colored_echo "$RED" "Workdir not defined in config file" >&2
	  colored_echo "$YELLOW" "Provide full path of a working directory: " >&2
	  read workdir
	  echo $workdir
  else
	  value="${value#\"}"
	  echo "${value%\"}"
  fi
}
