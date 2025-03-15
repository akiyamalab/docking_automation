===========
ドメイン層API
===========

ドメイン層は、システムの中核となるビジネスロジックを実装しています。
このセクションでは、ドメイン層の主要なモジュールとクラスのAPIリファレンスを提供します。

.. toctree::
   :maxdepth: 2

   molecule
   docking
   job
   result
   storage

分子準備ドメイン
==============

.. automodule:: docking_automation.domain.molecule
   :members:
   :undoc-members:
   :show-inheritance:

エンティティ
----------

.. automodule:: docking_automation.domain.molecule.entity
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: docking_automation.domain.molecule.entity.compound.Compound
   :members:
   :undoc-members:
   :show-inheritance:

値オブジェクト
-----------

.. automodule:: docking_automation.domain.molecule.value_object
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: docking_automation.domain.molecule.value_object.molecule_structure.MoleculeStructure
   :members:
   :undoc-members:
   :show-inheritance:

サービス
-------

.. automodule:: docking_automation.domain.molecule.service
   :members:
   :undoc-members:
   :show-inheritance:

リポジトリ
--------

.. automodule:: docking_automation.domain.molecule.repository
   :members:
   :undoc-members:
   :show-inheritance:

ドッキング実行ドメイン
==================

.. automodule:: docking_automation.domain.docking
   :members:
   :undoc-members:
   :show-inheritance:

エンティティ
----------

.. automodule:: docking_automation.domain.docking.entity
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: docking_automation.domain.docking.entity.receptor.Receptor
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: docking_automation.domain.docking.entity.ligand.Ligand
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: docking_automation.domain.docking.entity.docking_task.DockingTask
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: docking_automation.domain.docking.entity.docking_result.DockingResult
   :members:
   :undoc-members:
   :show-inheritance:

値オブジェクト
-----------

.. automodule:: docking_automation.domain.docking.value_object
   :members:
   :undoc-members:
   :show-inheritance:

サービス
-------

.. automodule:: docking_automation.domain.docking.service
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: docking_automation.domain.docking.service.docking_service.DockingService
   :members:
   :undoc-members:
   :show-inheritance:

ジョブ管理ドメイン
===============

.. automodule:: docking_automation.domain.job
   :members:
   :undoc-members:
   :show-inheritance:

結果管理ドメイン
============

.. automodule:: docking_automation.domain.result
   :members:
   :undoc-members:
   :show-inheritance:

データアクセスドメイン
==================

.. automodule:: docking_automation.domain.storage
   :members:
   :undoc-members:
   :show-inheritance: