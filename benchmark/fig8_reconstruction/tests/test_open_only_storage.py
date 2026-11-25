"""Regression tests for the user-requested v4 nonduplicated open-channel store."""
from __future__ import annotations
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
import h5py
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from primitive_store import (write_config_first, write_physics,load_primitives,
                             radial_grid,get_channel,get_forcing)


def build_test_record(tmp_path):
    c=json.loads((Path(__file__).resolve().parents[1]/'configs'/'fig8_point.json').read_text())
    c.update(log_points=5,orbit_interval_points=0)
    grid=radial_grid(c)
    record=lambda val:SimpleNamespace(coefficient=np.full(len(grid),val,dtype=complex),
                           quadrature_converged=np.ones(len(grid),dtype=bool))
    cache={(1,-1,1):record(0.5+1j),
           (1,-1,-1):record(123.0), # closed channel: must not be stored
           (3,2,1):record(999.0)}  # not used by open source: must not be stored
    source={(0,0,1):np.arange(len(grid),dtype=complex)+1j}
    amplitude={(0,0,1):.0001+.0007j}
    partial=[{'ell_f':0,'m_f':0,'g':1,'forcing_l':1,'I_l_re':0.25,'I_l_im':-0.05}]
    folder=tmp_path/'point'
    write_config_first(folder,c)
    write_physics(folder/'physics.h5',config=c,forcing_cache=cache,
                  source_cache=source,amplitude_cache=amplitude,
                  sampled_zero_keys=[],partial_records=partial)
    return folder,c,grid


def test_open_only_and_direct_indexing(tmp_path):
    folder,c,grid=build_test_record(tmp_path)
    assert sorted(p.name for p in folder.iterdir()) == ['config.json','physics.h5']
    assert json.loads((folder/'config.json').read_text())==c
    with h5py.File(folder/'physics.h5','r') as h:
        assert 'radius' not in h
        assert 'forcing/abs_error' not in h
        assert 'solver' not in h
        assert 'forcing/u_lmg' in h
        assert 'channels/S_total' in h and 'channels/Z_infinity' in h
        assert np.array_equal(h['forcing/keys_lmg'][:],[[1,-1,1]])
    d=load_primitives(folder)
    r,u=get_forcing(d,1,-1,1)
    assert np.array_equal(r,grid)
    assert np.all(u==.5+1j)
    with pytest.raises(KeyError):get_forcing(d,1,-1,-1)
    ch=get_channel(d,0,0,1)
    assert np.array_equal(ch['S'],np.arange(len(grid))+1j)
    assert ch['Z_infinity']==.0001+.0007j
    assert ch['I_l']=={1:.25-.05j}


def test_grid_checksum_rejects_silent_recipe_drift(tmp_path):
    folder,c,grid=build_test_record(tmp_path)
    path=folder/'config.json'
    config=json.loads(path.read_text())
    config['r_max']=2*config['r_max']
    path.write_text(json.dumps(config))
    with pytest.raises(ValueError,match='Recreated radial grid differs'):
        load_primitives(folder)


def test_existing_folder_rejected(tmp_path):
    folder,c,_=build_test_record(tmp_path)
    with pytest.raises(FileExistsError):
        write_config_first(folder,c)


def test_config_registered_before_physics_import(tmp_path):
    from one_point_fast import run
    c=json.loads((Path(__file__).resolve().parents[1]/'configs'/'fig8_point.json').read_text())
    dest=tmp_path/'failed_input_well_preserved'
    with pytest.raises(FileNotFoundError):
        run(c,tmp_path/'nonexistent_modules',dest)
    assert json.loads((dest/'config.json').read_text())==c
    assert not (dest/'physics.h5').exists()
