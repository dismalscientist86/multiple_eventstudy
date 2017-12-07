program define instanttex


syntax using/, [Table(string asis) LTable(string asis) Figure(string asis) force replace Package(string) NOcompile]




if `"`table'"' == `""' & `"`figure'"' == `""' & `"`ltable'"' == `""'{
    disp as err "Must Specific Either table(), ltable() or figure()"
    exit
}

if `"`table'"' != ""{
    foreach tab in `table'{
        if strpos("`tab'",".")==0 local tab = "`tab'" + ".tex"
        capture confirm file "`tab'"
        if _rc != 0{
            disp as err "File `tab' does not exist."
            if "`force'"!="force" exit
        }
    }
}

if `"`ltable'"' != ""{
    foreach tab in `ltable'{
        if strpos("`tab'",".")==0 local tab = "`tab'" + ".tex"
        capture confirm file "`tab'"
        if _rc != 0{
            disp as err "File `tab' does not exist."
            if "`force'"!="force" exit
        }
    }
}

if `"`figure'"' != "" {
    foreach fig in `figure'{
        if strpos("`fig'",".")==0 local fig = "`fig'" + ".eps"
        capture confirm file "`fig'"
        if _rc != 0{
            disp as err "File `fig' does not exist."
            if "`force'"!="force" exit
        }
    }
}


tempname output
file open `output' using "`using'", write `replace'
file write `output' "\documentclass{article}" _n "\usepackage{geometry, booktabs}" _n "\usepackage[pdftex]{graphicx}"_n "\usepackage[suffix=]{epstopdf}" _n "\usepackage{rotating}" _n

if "`package'" != "" {
    local first = 1
    file write `output' "\usepackage{"
    foreach pack in `package'{
        if `first' == 1 {
            file write `output' "`pack'"
            local `first' = 0
        }
        else file write `output' ", `pack'"
    }
    file write `output' "}" _n
}
if `"`ltable'"' != `""'{
    file write `output' "\usepackage{lscape}" _n
}

file write `output' "\begin{document}" _n 

if `"`table'"' !=""{
    foreach tab in `table'{
        if regexm("`tab'","[.]")==0 local tab = "`tab'" + ".tex"
        file write `output' "\input{`tab'} \pagebreak" _n
    }
}

if `"`ltable'"' !=""{
    foreach tab in `ltable'{
        if regexm("`tab'","[.]")==0 local tab = "`tab'" + ".tex"
        file write `output' "\begin{landscape}" _n "\input{`tab'}" _n "\end{landscape}" _n "\pagebreak" _n
    }
}

if `"`figure'"' != "" {
    foreach fig in `figure'{
        if regexm("`fig'","[.]")==0 local fig = "`fig'" + ".eps"
        file write `output' "\includegraphics[width=\textwidth]{`fig'} \pagebreak" _n
    }
}

file write `output' "\end{document}"
file close `output'

*local filename = regexr("`using'","^using","")

if "`nocompile'"==""{
    sh pdflatex -shell-escape "`using'"
}


end

