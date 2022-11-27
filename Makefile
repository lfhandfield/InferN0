
installed:
  if sh install.sh; then \
    echo success > installed; \
  else \
    echo 'failed to install, are you root?'; \
  fi
    
