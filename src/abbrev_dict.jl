#!/usr/bin/env julia

# println("\nDANS abbrevs_dict.jl")

# Inspiré du gem abbrev de Ruby
#   http://ruby-doc.org/stdlib-2.5.0/libdoc/abbrev/rdoc/Abbrev.html
#   https://github.com/ruby/ruby/blob/trunk/lib/abbrev.rb
#
# Given a set of strings, calculate the set of unambiguous abbreviations for
# those strings, and return a hash where the keys are all the possible
# abbreviations and the values are the full strings.
#
# Thus, given +words+ is "car" and "cone", the keys pointing to "car" would
# be "ca" and "car", while those pointing to "cone" would be "co", "con", and
# "cone".
#
#   abbrevsdict(["car", "cone"])
#   #=> {"ca"=>"car", "con"=>"cone", "co"=>"cone", "car"=>"car", "cone"=>"cone"}
#
# The optional +pattern+ parameter is a pattern or a string. Only input
# strings that match the pattern or start with the string are included in the
# output hash.
#
#   abbrevsdict(%w{car box cone crab}, /b/)
#   #=> {"box"=>"box", "bo"=>"box", "b"=>"box", "crab" => "crab"}
#
#   abbrevsdict(%w{car box cone}, 'ca')
#   #=> {"car"=>"car", "ca"=>"car"}
#
# TODO: en faire un Pkg Julia officiel
# TODO: accepter une collection de String au sens large
#       (e.g une Base.KeyIterator)
# TODO: créer un module AbbrevsDict et methode AbbrevsDict.abbrevsdict()
#
# HIST
# - xx/06/2018 création par améliorer ArgParse
# - 15/10/2018 table les clés deviennent des Strings plutot que des Symbols
function abbrevsdict(words::Array{String,1}, pat = r".*")
    if isa(pat, String)
        pat = Regex(pat)
    end
    # table = Dict{String, Symbol}()
    table = Dict{String,String}()
    seen = Dict{String,Int}()
    for word in words
        length(word) == 0 && next
        for len = length(word):-1:1
            abbrev = word[1:len]
            # if !ismatch(pat, abbrev)
            #     continue
            # end
            if !occursin(pat, abbrev)
                continue
            end
            # Without Base.xxx, could raise error if get is redefined somewhere
            seen[abbrev] = 1 + Base.get(seen, abbrev, 0)
            if seen[abbrev] == 1
                # table[abbrev] = Symbol(word)
                table[abbrev] = word
            elseif seen[abbrev] == 2
                delete!(table, abbrev)
            else
                break
            end
        end
    end
    for word in words
        if !occursin(pat, word)
            continue
        end
        # table[word] = Symbol(word)
        table[word] = word
    end
    table
end
# function abbrevsdict(words::Array{String,1}, pat=r".*")
#     if isa(pat, String)
#         pat = Regex(pat)
#     end
#     table = Dict{String, Symbol}()
#     seen = Dict{String, Int}()
#     for word in words
#         length(word)==0 && next
#         for len in length(word):-1:1
#             abbrev = word[1:len]
#             @show abbrev, typeof(abbrev)
#             if !occursin(pat, abbrev)
#                 continue
#             end
#             # Without Base.xxx, could raise error if get is redefined somewhere
#             seen[abbrev] = 1 + Base.get(seen, abbrev, 0)
#             if seen[abbrev] == 1
#                 table[abbrev] = Symbol(word)
#             elseif seen[abbrev] == 2
#                 delete!(table, abbrev)
#             else
#                 break
#             end
#         end
#     end
#     @show words
#     # for word in words
#     #     # if !ismatch(pat, word)
#     #     #     continue
#     #     # end
#     #     if !occursin(pat, word)
#     #         continue
#     #     end
#     #     table[word] = word
#     # end
#     @show table
#     table
# end

function test_abbrevsdict()
    @show typeof(ARGS)
    @show ARGS
    @show abbrevsdict(ARGS)
    @show abbrevsdict(["car", "cone"])
    @show abbrevsdict(["car", "cone"], "e")
    @show abbrevsdict(["car", "cone"], r"[q-z]")
end

if basename(@__FILE__) == basename(PROGRAM_FILE)
    test_abbrevsdict()
    # main()
end
