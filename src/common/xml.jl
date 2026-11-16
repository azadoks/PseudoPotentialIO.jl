"""
Some UPF files seem to have unescaped ampersands which cause XML parsing to fail.

This could also be a problem for other XML-based formats (e.g. PSML).

In XML, `&` actually _begins_ escape sequences (e.g. `<` is `&lt;`),
so we only want to replace `&`s that are not followed by a valid escape sequence.
"""
escape_ampersand(text) = replace(text, r"&(?!(amp|apos|quot|lt|gt);)" => "&amp;")

"""
Print the XML of `node` (and its descendants), indented, to `io`: the output of
`EzXML.prettyprint(io, node)`.

Thread-safe, unlike `EzXML.prettyprint` (EzXML 1.2.3): that one dumps into a libxml2
buffer owned by a Julia object with a finalizer and does not keep the object alive while
it copies the text out, so a garbage collection started by any thread at that moment
(the copy allocates the `String`) frees the buffer under it, and the text read is freed
memory (PseudoPotentialIO #7). Here the buffer is freed explicitly, after the copy, and
`node` is kept alive throughout.
"""
function xml_prettyprint(io::IO, node::EzXML.Node)
    doc_ptr = EzXML.hasdocument(node) ? EzXML.document(node).node.ptr : C_NULL
    GC.@preserve node begin
        buf = ccall((:xmlBufferCreate, EzXML.libxml2), Ptr{Cvoid}, ())
        buf == C_NULL && throw(OutOfMemoryError())
        try
            # level 0, format 1: indented, as `EzXML.prettyprint`
            len = ccall((:xmlNodeDump, EzXML.libxml2), Cint,
                        (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Cint, Cint),
                        buf, doc_ptr, node.ptr, 0, 1)
            len == -1 && error("libxml2 failed to dump the XML node $(nodename(node))")
            content = ccall((:xmlBufferContent, EzXML.libxml2), Ptr{UInt8}, (Ptr{Cvoid},),
                            buf)
            nbytes = ccall((:xmlBufferLength, EzXML.libxml2), Cint, (Ptr{Cvoid},), buf)
            print(io, unsafe_string(content, nbytes))
        finally
            ccall((:xmlBufferFree, EzXML.libxml2), Cvoid, (Ptr{Cvoid},), buf)
        end
    end
    return nothing
end
