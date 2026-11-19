# Source provenance for oracle reports: which commits of Jutul, JutulDarcy and
# MRST a comparison ran against.
#
# Why: Jutul and JutulDarcy may be `dev`-ed from local checkouts, for which
# Manifest.toml records no commit. So every `OracleReport` records the commit
# and dirty state of each checkout when it is made, or the version of a
# registered release.
#
# Record-only by design: nothing here fails a gate; a dirty checkout gives a
# warning, not an error.

# -------------------------------- Source trees --------------------------------

"""
    SourceInfo

Provenance of one source tree (Jutul, JutulDarcy or MRST), as recorded in
an `OracleReport`. Shown as `"name commit"`, plus `" (dirty)"` if dirty.

# Fields

- `name`: the source tree's name (`"Jutul"`, `"JutulDarcy"` or `"MRST"`).
- `dir`: the checkout directory; empty when unknown or not loaded.
- `commit`: the short (8-character) commit hash of `HEAD` for a git
  checkout; `"v<version>"` for a registered release of a loaded package;
  `"unknown"` when neither is known, `"not loaded"` for a package that is not
  loaded.
- `dirty`: whether tracked files have uncommitted changes; `missing` when
  unknown.
"""
struct SourceInfo
    name::String                 # package
    dir::String                  # checkout directory
    commit::String               # short commit hash, or v<version> of a release
    dirty::Union{Bool, Missing}  # uncommitted changes (missing if unknown)
end

function Base.show(io::IO, s::SourceInfo)
    print(io, s.name, " ", s.commit)
    s.dirty === true && print(io, " (dirty)")
end

# Packages whose checkouts are tracked, looked up among the loaded modules by
# UUID so that Oracle does not need to depend on them.
const TRACKED_PACKAGES = (
    Base.PkgId(Base.UUID("2b460a1a-8a2b-45b2-b125-b5c536396eb9"), "Jutul"),
    Base.PkgId(Base.UUID("82210473-ab04-4dce-b31b-11573c4f8e0a"), "JutulDarcy"),
)

# Output of `git -C dir args...`, stripped; `nothing` if git fails or is missing
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

# Returns

A new vector of [`SourceInfo`](@ref), Jutul first. A git checkout has its
commit; a registered release its version (`"v0.5.9"`), not dirty. A package
that is not loaded has commit `"not loaded"` and `dirty` `missing`.
"""
function source_provenance()
    infos = SourceInfo[]
    for id in TRACKED_PACKAGES
        if haskey(Base.loaded_modules, id)
            mod = Base.loaded_modules[id]
            info = _git_source_info(id.name, pkgdir(mod))
            # not a git checkout: a registered release, identified by its version
            if info.commit == "unknown" && !isnothing(pkgversion(mod))
                info = SourceInfo(id.name, info.dir, "v$(pkgversion(mod))", false)
            end
            push!(infos, info)
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

# Warn for each source tree whose checkout has uncommitted changes
function _warn_dirty(infos)
    for s in infos
        s.dirty === true &&
            @warn "$(s.name) checkout has uncommitted changes to tracked files; " *
                  "this report cannot be reproduced from commit $(s.commit) alone" dir =
                  s.dir
    end
end
