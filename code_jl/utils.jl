function TexTable(name::String, headings1::Vector{String}, headings2::Vector{String}, body::Array{T, 2} where T <: Real, types::Vector{DataType}, digits = 3, head_cnt = 2)

    @assert size(headings1, 1) == size(headings2, 1)
    @assert size(headings1, 1) == size(body, 2)

    @inline function _write_heading(head::Vector{String})

        len = size(head, 1)
        last_seen = "nAn"
        dub_counter = 1
        counter = 1
        for h in head
            if counter < len
                if h != last_seen && head[counter + 1] != h
                    if dub_counter == 1
                        write(io, h)
                        write(io, "&")
                        last_seen = h
                    else
                        write(io, "\\multicolumn{$dub_counter}{c}{$last_seen}")
                        write(io, "&")
                        last_seen = h
                        dub_counter = 1
                    end
                else
                    dub_counter += 1
                end
            else
                if dub_counter != 1 || h == last_seen
                    write(io, "\\multicolumn{$dub_counter}{c}{$last_seen}")
                    write(io, "\\\\")
                    last_seen = h
                else
                    write(io, h)
                    write(io, "\\\\")
                    last_seen = h
                end
            end
            last_seen = h
            counter += 1
        end
    end

    @inline function _write_body(body::Array{T, 2} where T <: Any)

        rows = size(body, 1)
        cols = size(body, 2)

        for r in 1:rows
            for c in 1:cols
                if types[c] == Float64
                    write(io, string(round(abs(body[r, c]), digits = digits)))
                elseif types[c] == Int64
                    write(io, string(convert(Int64, body[r, c])))
                end
                if c == cols
                    write(io, "\\\\\n")
                else
                    write(io, "&")
                end
            end
        end
    end

    align = ""
    for i in headings1
        align = string(align, "l")
    end

    io = open(name, "w")
    write(io, "\\begin{table}\n")
    write(io, "\\begin{tabular}{$align}\n")
    write(io, "\\toprule\n")
    _write_heading(headings1)
    write(io, "\n\\midrule\n")
    if head_cnt == 2
        _write_heading(headings2)
        write(io, "\n\\toprule\n")
    end
    _write_body(body)
    write(io, "\\bottomrule\n")
    write(io, "\\end{tabular}\n")
    write(io, "\\end{table}\n")
    close(io)

end
