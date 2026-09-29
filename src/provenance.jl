# Source provenance for oracle reports.
#
# Jutul and JutulDarcy are `dev`-ed from local checkouts, and for a dev-ed
# package Manifest.toml records only `path` and `version` — no commit, no tree
# hash. It therefore does not pin what a gate actually ran against. This file
# fills that gap: every `OracleReport` captures the commit and dirty state of
# each checkout at construction time, so a result can always be attributed to
# the source trees that produced it.
#
# Record-only by design: nothing here fails a gate. A dirty checkout triggers a
# warning, not an error.

# -------------------------------- Source trees --------------------------------

"""
    SourceInfo

Provenance of one source tree: package `name`, checkout `dir`, short `commit`
hash, and whether tracked files have uncommitted changes (`dirty`).

`commit` is `"unknown"` when the directory is not a git repository or git is
unavailable; `dirty` is then `missing`.
"""
struct SourceInfo
    name::String                 # package
    dir::String                  # checkout directory
    commit::String               # short commit hash ("unknown" if none)
    dirty::Union{Bool, Missing}  # uncommitted changes (missing if unknown)
end

function Base.show(io::IO, s::SourceInfo)
    print(io, s.name, " ", s.commit)
    s.dirty === true && print(io, " (dirty)")
end

# Packages whose checkouts are tracked. UUIDs from their Project.toml; looked up
# among loaded modules so that Oracle does not need to depend on them.
const TRACKED_PACKAGES = (
    Base.PkgId(Base.UUID("2b460a1a-8a2b-45b2-b125-b5c536396eb9"), "Jutul"),
    Base.PkgId(Base.UUID("82210473-ab04-4dce-b31b-11573c4f8e0a"), "JutulDarcy"),
)

function _git(dir::AbstractString, args::Vector{String})
    cmd = Cmd(`git -C $dir $args`; ignorestatus = true)
    out = IOBuffer()
    p = try
        run(pipeline(cmd; stdout = out, stderr = devnull))
    catch
        return nothing          # git not installed
    end
    success(p) || return nothing
    return strip(String(take!(out)))
end

"""
    _git_source_info(name, dir) -> SourceInfo

Commit and dirty state of the checkout at `dir`.

Dirty means uncommitted changes to *tracked* files. Untracked files are
ignored: they cannot change loaded code unless a tracked file includes them,
in which case that tracked file shows as modified.
"""
function _git_source_info(name::AbstractString, dir::AbstractString)
    commit = _git(dir, ["rev-parse", "--short=8", "HEAD"])
    commit === nothing && return SourceInfo(name, dir, "unknown", missing)

    status = _git(dir, ["status", "--porcelain", "--untracked-files=no"])
    dirty = status === nothing ? missing : !isempty(status)
    return SourceInfo(name, dir, commit, dirty)
end

"""
    source_provenance() -> Vector{SourceInfo}

Provenance of the Jutul and JutulDarcy source trees currently loaded.

A package that is not loaded is reported with commit `"not loaded"`.
"""
function source_provenance()
    infos = SourceInfo[]
    for id in TRACKED_PACKAGES
        if haskey(Base.loaded_modules, id)
            dir = pkgdir(Base.loaded_modules[id])
            push!(infos, _git_source_info(id.name, dir))
        else
            push!(infos, SourceInfo(id.name, "", "not loaded", missing))
        end
    end
    return infos
end

# ------------------------- MRST, and dirty checkouts --------------------------

"""
    _mrst_source_info(dump) -> SourceInfo

MRST provenance as recorded by the MATLAB side (`mrst_commit`, `mrst_root`).
Dumps without those fields report `"unknown"`.
"""
function _mrst_source_info(dump)
    commit = haskey(dump, "mrst_commit") ? String(dump["mrst_commit"]) : "unknown"
    root   = haskey(dump, "mrst_root")   ? String(dump["mrst_root"])   : ""
    commit != "unknown" && (commit = first(commit, 8))
    return SourceInfo("MRST", root, commit, missing)
end

function _warn_dirty(infos)
    for s in infos
        s.dirty === true &&
            @warn "$(s.name) checkout has uncommitted changes to tracked files; " *
                  "this report cannot be reproduced from commit $(s.commit) alone" dir =
                  s.dir
    end
end
