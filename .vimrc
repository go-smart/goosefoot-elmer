setlocal tabstop=4
setlocal shiftwidth=4
setlocal expandtab
setlocal autoindent
setlocal smarttab
setlocal formatoptions=croql

autocmd BufNewFile,BufRead *.src set filetype=fortran

if (@% =~ '\.src$')
  let fortran_free_source=1
  let fortran_more_precise=1
  set filetype=fortran
endif

let g:syntastic_python_flake8_args='--ignore=E501'
